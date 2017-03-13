/* sdsl - succinct data structures library
 *    Copyright (C) 2009 Simon Gog
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see http://www.gnu.org/licenses/ .
 */
/*! \file rmq_succinct_rec_new.hpp
 *    \brief rmq_succinct_rec_new.hpp contains the class rmq_succinct_bp_fast which supports range minimum or range maximum queries on a random access container in constant time and \f$2 n+o(n) bits\f$ space.
 *    \author Tobias Heuer
 */
#ifndef INCLUDED_SDSL_RMQ_SUCCINCT_REC_NEW
#define INCLUDED_SDSL_RMQ_SUCCINCT_REC_NEW

#include <stack>
#include <limits>

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "bits.hpp"
#include "bp_support_sada.hpp"
#include "bp_support_algorithm.hpp"
#include "rank_select_support_bp.hpp"
#include "suffix_tree_helper.hpp"
#include "util.hpp"


//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint32_t t_super_block_size=1024,uint32_t t_block_size=128, uint32_t t_levels=2, bool t_min = true>
class rmq_succinct_rec_new;

template<uint32_t t_super_block_size=1024, uint32_t t_block_size=128, uint32_t t_levels=2, bool t_min = false>
struct range_maximum_bp_fast_rec_new {
    typedef rmq_succinct_rec_new<t_super_block_size, t_block_size, t_levels, t_min> type;
};

//! A class to support range minimum or range maximum queries on a random access container.
/*!
 *  \tparam t_min        Specifies whether the data structure should answer range min/max queries (mimumum=true)
 *  \tparam t_bp_support Type of Support structure for the BPS-SCT.
 *
 * \par Time complexity
 *        \f$ \Order{1} \f$ for the range minimum/maximum queries if the balanced parentheses support structure supports constant time operations.
 * \par Space complexity:
 *        \f$ \Order{2n}+o(n) \f$ bits for the data structure ( \f$ n=size() \f$ ).
 *
 * \par Reference
 * H . Ferrada and G. Navarro.
 * Improved Range Minimum Queries.
 * In Proceedings of Data Compression Conference, DCC'16.
 *
 */
template<uint32_t t_super_block_size, uint32_t t_block_size, uint32_t t_levels, bool t_min>
class rmq_succinct_rec_new
{

        bool                        m_use_sparse_rmq;
        bit_vector                  m_gct_bp;         //!< A bit vector which contains the BP-GT of the input container.
        int_vector<>                m_min_excess_idx;
        int_vector<>                m_min_excess;
        rmq_support_sparse_table<true,false>  m_sparse_rmq;
        rmq_succinct_rec_new<t_block_size,t_block_size, (t_levels == 0 ? 0 : t_levels-1)>*  m_rmq_recursive;
        rank_select_support_bp<>      m_rank_select;
        bit_vector::value_type      m_max_excess_v;
        bit_vector::value_type      m_max_excess_reverse_v;

        void copy(const rmq_succinct_rec_new& rm) {
            m_use_sparse_rmq = rm.m_use_sparse_rmq;
            if (!m_use_sparse_rmq) {            
                m_gct_bp = rm.m_gct_bp;
                m_rank_select = rm.m_rank_select;
                m_rank_select.set_vector(&m_gct_bp);
                m_min_excess = rm.m_min_excess;
                m_min_excess_idx = rm.min_excess_idx;
                m_rmq_recursive = rm.m_rmq_recursive;
                m_max_excess_v = rm.m_max_excess_v;
                m_max_excess_reverse_v = rm.m_max_excess_reverse_v;
            } else {
                m_min_excess = rm.m_min_excess;
                m_sparse_rmq = rm.m_sparse_rmq;
                m_sparse_rmq.set_vector(&m_min_excess);
            }
        }

    private:
        template<class t_rac>
        rmq_succinct_rec_new() : m_rmq_recursive(nullptr) { }


        template<bool t_strict, bool t_reverse, class t_rac>
        bit_vector::value_type construct_generalized_cartesian_tree(const t_rac* v, bool write_bp_sequence=true) {
            typedef min_max_trait<t_rac, true, t_strict> mm_trait;
            bit_vector::value_type max_excess = 0, cur_excess = 0;
            if (v->size() > 0) {
                long int cur_pos = t_reverse ? v->size()-1 : 0;
                size_t bp_cur_pos = 0;
                std::stack<typename t_rac::value_type> s;
                s.push(std::numeric_limits<typename t_rac::value_type>::min());
                if (write_bp_sequence) m_gct_bp[bp_cur_pos++] = 1;
                while ((!t_reverse && cur_pos <  ((long int) v->size())) || (t_reverse && cur_pos >= 0)) {
                    typename t_rac::value_type cur_elem = (*v)[cur_pos];
                    cur_pos += t_reverse ? -1 : 1;
                    while (mm_trait::compare(cur_elem, s.top()) && s.size() > 1) {
                        s.pop();
                        bp_cur_pos++; cur_excess--;
                    }
                    if (write_bp_sequence) m_gct_bp[bp_cur_pos++] = 1;
                    cur_excess++;
                    if (cur_excess > max_excess) max_excess = cur_excess;
                    s.push(cur_elem);
                }
                while (!s.empty()) {
                    s.pop();
                    bp_cur_pos++;
                }
            }
            return max_excess;
        }

        void build_rmq_recursive() {
            size_type bp_size = m_gct_bp.size();
            size_type excess_block_size = bp_size/t_super_block_size + (bp_size % t_super_block_size != 0 ? 1 : 0); 
            m_min_excess = int_vector<>(excess_block_size,0);
            m_min_excess_idx = int_vector<>(excess_block_size,0);
            bit_vector::difference_type min_rel_ex = 0;
            for (size_t i = 0; i*t_super_block_size < bp_size; ++i) {
                uint64_t min_idx = near_rmq(m_gct_bp,i*t_super_block_size, std::min((i+1)*t_super_block_size - 1,bp_size-1),min_rel_ex);
                m_min_excess_idx[i] = min_idx-i*t_super_block_size;
                m_min_excess[excess_block_size - 1 - i] = m_rank_select.excess(min_idx);
            }
            m_rmq_recursive = new rmq_succinct_rec_new<t_block_size,t_block_size,(t_levels == 0 ? 0 : t_levels-1)>(&m_min_excess);
            m_use_sparse_rmq = false;
            util::bit_compress(m_min_excess);
            util::bit_compress(m_min_excess_idx);
        }

        inline std::pair<bit_vector::size_type, int_vector<>::value_type> min_ex(const bit_vector::size_type i1, const bit_vector::size_type i2, const bit_vector::size_type i3,
                const int_vector<>::value_type i1_ex, const int_vector<>::value_type i2_ex, const int_vector<>::value_type i3_ex) const {
            assert(i1 <= i2); assert(i2 <= i3);
            if (i3_ex <= i1_ex && i3_ex <= i2_ex) {
                return std::make_pair(i3,i3_ex);
            } else if(i2_ex <= i1_ex) {
                return std::make_pair(i2,i2_ex);
            } else {
                return std::make_pair(i1,i1_ex);
            } 
        }

        inline int_vector<>::value_type get_min_excess_idx(size_t i) const {
            return m_min_excess_idx[i] + i * t_super_block_size;
        }

        inline int_vector<>::size_type map_to_min_excess(size_t i) const {
            return m_min_excess.size()-1-i;
        }

        inline int_vector<>::size_type get_min_excess(size_t i) const {
            return m_min_excess[map_to_min_excess(i)];
        }

        inline int_vector<>::size_type map_index(size_t i) {
            if(m_max_excess_v < m_max_excess_reverse_v) {
                return i;
            } else {
                return size() - 1 - i;
            }
        }

    public:
        typedef typename bit_vector::size_type size_type;
        typedef typename bit_vector::value_type value_type;
        typedef typename int_vector<>::value_type i_value_type;

        //! Default constructor
        rmq_succinct_rec_new() {}

        //! Constructor
        /*! \tparam t_rac A random access container.
         *  \param  v     Pointer to container object.
         */
        template<class t_rac>
        rmq_succinct_rec_new(const t_rac* v=nullptr) : m_rmq_recursive(nullptr) {
            if (v != nullptr) {
                size_type bp_size = 2*v->size()+2;
                m_gct_bp = bit_vector(bp_size,0);
                if(t_levels > 0 && t_super_block_size <= bp_size) {
                    m_max_excess_v = construct_generalized_cartesian_tree<true,false>(v,false);
                    m_max_excess_reverse_v = construct_generalized_cartesian_tree<false,true>(v,false);
                    if (m_max_excess_v < m_max_excess_reverse_v) {
                        construct_generalized_cartesian_tree<true,false>(v);
                    } else {
                        construct_generalized_cartesian_tree<false,true>(v);
                    }
                    m_rank_select = rank_select_support_bp<>(&m_gct_bp);
                    build_rmq_recursive();
                }
                else {
                    m_min_excess = int_vector<>(v->size(),0);
                    for (size_t i = 0; i < v->size(); ++i)
                    {
                        m_min_excess[i] = (*v)[i];
                    }
                    m_sparse_rmq = rmq_support_sparse_table<true,false>(&m_min_excess);
                    m_use_sparse_rmq = true;
                }
            }
        }

        ~rmq_succinct_rec_new() {
            if (t_levels > 0 && m_gct_bp.size() > 2) {
                delete m_rmq_recursive;
            }
        }

        //! Copy constructor
        rmq_succinct_rec_new(const rmq_succinct_rec_new& rm) {
            *this = rm;
        }

        //! Move constructor
        rmq_succinct_rec_new(rmq_succinct_rec_new&& rm) {
            *this = std::move(rm);
        }

        rmq_succinct_rec_new& operator=(const rmq_succinct_rec_new& rm) {
            if (this != &rm) {
                m_use_sparse_rmq = rm.m_use_sparse_rmq;
                if (!m_use_sparse_rmq) {
                    m_gct_bp = rm.m_gct_bp; 
                    m_rank_select = rm.m_rank_select;
                    m_rank_select.set_vector(&m_gct_bp);
                    m_min_excess = rm.m_min_excess;
                    m_min_excess_idx = rm.min_excess_idx;
                    m_rmq_recursive = rm.m_rmq_recursive;
                    m_max_excess_v = rm.m_max_excess_v;
                    m_max_excess_reverse_v = rm.m_max_excess_reverse_v;
                } else {
                    m_min_excess = rm.m_min_excess;
                    m_sparse_rmq = rm.m_sparse_rmq;
                    m_sparse_rmq.set_vector(&m_min_excess);
                }
            }
            return *this;
        }

        rmq_succinct_rec_new& operator=(rmq_succinct_rec_new&& rm) {
            if (this != &rm) {
                m_use_sparse_rmq = rm.m_use_sparse_rmq;
                if (!m_use_sparse_rmq) {
                    m_gct_bp = std::move(rm.m_gct_bp);
                    m_rank_select = std::move(rm.m_rank_select);
                    m_rank_select.set_vector(&m_gct_bp);
                    m_min_excess = std::move(rm.m_min_excess);
                    m_min_excess_idx = std::move(rm.m_min_excess_idx);
                    m_rmq_recursive = std::move(rm.m_rmq_recursive);
                    m_max_excess_v = rm.m_max_excess_v;
                    m_max_excess_reverse_v = rm.m_max_excess_reverse_v;
                } else {
                    m_min_excess = std::move(rm.m_min_excess);
                    m_sparse_rmq = std::move(rm.m_sparse_rmq);
                    m_sparse_rmq.set_vector(&m_min_excess);
                }
            }
            return *this;
        }

        void swap(rmq_succinct_rec_new& rm) {
            std::swap(m_use_sparse_rmq, rm.m_use_sparse_rmq);
            if (!m_use_sparse_rmq) {
                m_gct_bp.swap(rm.m_gct_bp);
                util::swap_support(m_rank_select, rm.m_rank_select,
                                &m_gct_bp, &(rm.m_gct_bp));
                m_min_excess.swap(rm.m_min_excess);
                m_min_excess_idx.swap(rm.m_min_excess_idx);
                *m_rmq_recursive.swap(*rm.m_rmq_recursive);
                std::swap(m_max_excess_v, rm.m_max_excess_v);
                std::swap(m_max_excess_reverse_v, rm.m_max_excess_reverse_v);
            } else {
                m_min_excess.swap(rm.m_min_excess);
                util::swap_support(m_sparse_rmq, rm.m_sparse_rmq,
                                   &m_min_excess, &(rm.min_excess));
            }
        }


        //! Range minimum/maximum query for the supported random access container v.
        /*!
         * \param l Leftmost position of the interval \f$[\ell..r]\f$.
         * \param r Rightmost position of the interval \f$[\ell..r]\f$.
         * \return The minimal index i with \f$\ell \leq i \leq r\f$ for which \f$ v[i] \f$ is minimal/maximal.
         * \pre
         *   - r < size()
         *   - \f$ \ell \leq r \f$
         * \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type operator()(const size_type l, const size_type r) {
            assert(l <= r); assert(r < size());
            if (l == r) return l;
            if(m_use_sparse_rmq) {
                 return m_sparse_rmq(l,r);
            }

            size_type tmp_l = map_index(l), tmp_r = map_index(r);
            if (tmp_l > tmp_r) std::swap(tmp_l,tmp_r);

            size_type i = m_rank_select.select(tmp_l+2)-1;

            /*if (tmp_r - tmp_l < 64) {
                value_type data = m_gct_bp.get_int(i+1);
                uint64_t one_cnt = bits::cnt(data)-1;
                if (tmp_l + one_cnt >= tmp_r) {
                    size_t cur_select = tmp_l+1, min_idx = i;
                    int cur_excess = ((tmp_l+1) << 1)-(i+1), min_excess = cur_excess;
                    for (size_t k = i+1; k <= i+65; ++k) {
                        if (m_gct_bp[k]) {
                            cur_excess++; cur_select++;
                            if (cur_select == tmp_r + 2) break;
                        } else {
                            cur_excess--;
                            if (cur_excess <= min_excess) {
                                min_excess = cur_excess;
                                min_idx = k;
                            }
                        }
                    }

                    if (t_block_size == 0 || m_max_excess_v < m_max_excess_reverse_v) {
                        return (min_excess+min_idx)>>1;
                    } else {
                        return N-((min_excess+min_idx)>>1)-1;
                    }
                }
            }*/
            size_type j = m_rank_select.select(tmp_r+2,i);

            size_type block_i = (i+t_super_block_size-1)/t_super_block_size;
            size_type block_j = j/t_super_block_size;
            
            bit_vector::difference_type min_rel_ex = 0;
            if (block_i >= block_j) {
                block_i = block_i-(i != 0);
                size_type min_excess_idx = i;
                if (block_i == block_j) {
                    min_excess_idx = get_min_excess_idx(block_i);
                } else {
                    size_type min_left_excess_idx = get_min_excess_idx(block_i);
                    size_type min_right_excess_idx = get_min_excess_idx(block_j);
                    min_excess_idx = (get_min_excess(block_i) < get_min_excess(block_j) ? min_left_excess_idx : min_right_excess_idx);
                }
                if (min_excess_idx < i || min_excess_idx > j) {
                    min_excess_idx = near_rmq(m_gct_bp,i,j-1,min_rel_ex); 
                }
                return map_index(m_rank_select.rank(min_excess_idx+1)-1);
            } else {
                size_type min_block_idx = (*m_rmq_recursive)(map_to_min_excess(block_j-1), map_to_min_excess(block_i));
                size_type min_block_excess_idx = get_min_excess_idx(map_to_min_excess(min_block_idx));
                i_value_type min_block_excess = m_min_excess[min_block_idx];
                
                size_type min_left_excess_idx = get_min_excess_idx(block_i-(i != 0)); 
                i_value_type min_left_excess = get_min_excess(block_i-(i != 0));
                if (min_left_excess < min_block_excess && min_left_excess_idx < i) {
                    min_left_excess_idx = near_rmq(m_gct_bp,i,t_super_block_size*block_i,min_rel_ex); 
                    min_left_excess = m_rank_select.excess(min_left_excess_idx); 
                }

                size_type min_right_excess_idx = get_min_excess_idx(block_j); 
                i_value_type min_right_excess = get_min_excess(block_j);
                if (min_right_excess <= min_block_excess && min_right_excess_idx > j) {
                    min_right_excess_idx = near_rmq(m_gct_bp,t_super_block_size*block_j,j,min_rel_ex); 
                    min_right_excess = m_rank_select.excess(min_right_excess_idx);
                }

                auto rmq_min = min_ex(min_left_excess_idx,min_block_excess_idx,min_right_excess_idx,
                                      min_left_excess,min_block_excess,min_right_excess);
                size_type min_idx = rmq_min.first;
                int_vector<>::value_type min_ex = rmq_min.second;
                return map_index((min_ex+min_idx)>>1);
            }
        }


        size_type size()const {
            if(m_use_sparse_rmq) {
                return m_sparse_rmq.size();
            }
            else {
              return (m_gct_bp.size()-2)/2;
            }
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_use_sparse_rmq, out, child, "m_use_sparse_rmq");
            if(m_rmq_recursive != nullptr) {
                written_bytes += m_gct_bp.serialize(out, child, "gct_bp");
                written_bytes += m_rank_select.serialize(out, child, "rank_select_bp");
                written_bytes += m_min_excess.serialize(out, child, "min_excess");
                written_bytes += m_min_excess_idx.serialize(out, child, "min_excess_idx");
                written_bytes += m_rmq_recursive->serialize(out, child, "rmq_recursive");
                written_bytes += write_member(m_max_excess_v, out, child, "m_max_excess_v");
                written_bytes += write_member(m_max_excess_reverse_v, out, child, "m_max_excess_v_reverse");
            } else {
                written_bytes += m_min_excess.serialize(out, child, "min_excess");
                written_bytes += m_sparse_rmq.serialize(out, child, "sparse_rmq");
            }
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            read_member(m_use_sparse_rmq,in);
            if (!m_use_sparse_rmq) {
                m_gct_bp.load(in);
                m_rank_select.load(in, &m_gct_bp);
                m_min_excess.load(in);
                m_min_excess_idx.load(in);
                m_rmq_recursive = new rmq_succinct_rec_new<t_block_size,t_block_size,(t_levels == 0 ? 0 : t_levels-1)>();
                m_rmq_recursive->load(in);
                read_member(m_max_excess_v,in);
                read_member(m_max_excess_reverse_v,in);
            } else {
                m_min_excess.load(in);
                m_sparse_rmq.load(in,&m_min_excess);
            }
        }
};

} // end namespace sdsl
#endif

