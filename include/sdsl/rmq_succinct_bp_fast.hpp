/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file rmq_succinct_bp_fast.hpp
    \brief rmq_succinct_bp_fast.hpp contains the class rmq_succinct_bp_fast which supports range minimum or range maximum queries on a random access container in constant time and \f$2 n+o(n) bits\f$ space.
    \author Tobias Heuer
*/
#ifndef INCLUDED_SDSL_RMQ_SUCCINCT_BP_FAST
#define INCLUDED_SDSL_RMQ_SUCCINCT_BP_FAST

#include <stack>
#include <limits>

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "bp_support_sada.hpp"
#include "suffix_tree_helper.hpp"
#include "util.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint32_t t_block_size=4096, bool t_min = true,
         class t_bp_support = bp_support_sada<t_block_size,32,rank_support_v5<> > >
class rmq_succinct_bp_fast;

template<uint32_t t_block_size=4096, class t_bp_support = bp_support_sada<t_block_size,32,rank_support_v5<> > >
struct range_maximum_bp_fast {
    typedef rmq_succinct_bp_fast<t_block_size, false, t_bp_support> type;
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
template<uint32_t t_block_size, bool t_min, class t_bp_support>
class rmq_succinct_bp_fast
{
        bit_vector                  m_gct_bp;         //!< A bit vector which contains the BP-GT of the input container.
        t_bp_support                m_gct_bp_support; //!< Support structure for the BPS-GT
        rmq_support_sparse_table<>  m_sparse_rmq;
        int_vector<>                m_min_excess_idx;
        int_vector<>                m_min_excess;        

        void copy(const rmq_succinct_bp_fast& rm) {
            m_gct_bp = rm.m_gct_bp;
            m_gct_bp_support = rm.m_gct_bp_support;
            m_gct_bp_support.set_vector(&m_gct_bp);
            m_sparse_rmq = rm.m_sparse_rmq;
        }

        
    private:
        
        template<class t_rac>
        void construct_generalized_cartesian_tree_rightmost(const t_rac* v) {
            if(v->size() > 0) {
                int64_t cur_pos = v->size()-1;
                int64_t bp_cur_pos = m_gct_bp.size()-1;
                std::stack<typename t_rac::value_type> s;
                s.push(std::numeric_limits<typename t_rac::value_type>::min()); 
                m_gct_bp[bp_cur_pos--] = 0;
                while(cur_pos >= 0) {
                    typename t_rac::value_type cur_elem = (*v)[cur_pos--];
                    while(s.top() >= cur_elem) {
                        s.pop();
                        m_gct_bp[bp_cur_pos--] = 1;
                    }
                    bp_cur_pos--;
                    s.push(cur_elem);
                }
                while(!s.empty()) {
                    s.pop();
                    m_gct_bp[bp_cur_pos--] = 1;
                }
            }
        }
        
        template<class t_rac>
        void construct_generalized_cartesian_tree_leftmost(const t_rac* v) {
            if(v->size() > 0) {
                size_t cur_pos = 0;
                size_t bp_cur_pos = 0;
                std::stack<typename t_rac::value_type> s;
                s.push(std::numeric_limits<typename t_rac::value_type>::min()); 
                m_gct_bp[bp_cur_pos++] = 1;
                while(cur_pos < v->size()) {
                    typename t_rac::value_type cur_elem = (*v)[cur_pos++];
                    while(s.top() > cur_elem) {
                        s.pop();
                        bp_cur_pos++;
                    }
                    m_gct_bp[bp_cur_pos++] = 1;
                    s.push(cur_elem);
                }
                while(!s.empty()) {
                    s.pop();
                    bp_cur_pos++;
                }
            }
        }
        
        void construct_sparse_table_over_sampled_excess_value_of_bp() {
            m_min_excess = int_vector<>(m_gct_bp.size()/t_block_size,0);
            m_min_excess_idx = int_vector<>(m_gct_bp.size()/t_block_size,0);
            for(size_t i = 0; (i+1)*t_block_size < m_gct_bp.size(); ++i) {
                uint64_t min_idx = m_gct_bp_support.rmq(i*t_block_size, (i+1)*t_block_size - 1);
                m_min_excess_idx[i] = min_idx;
                m_min_excess[i] = m_gct_bp_support.excess(min_idx);
            }
            m_sparse_rmq = rmq_support_sparse_table<>(&m_min_excess);
        }
        
        inline bit_vector::size_type min_excess(const bit_vector::size_type i1, const bit_vector::size_type i2, const bit_vector::size_type i3) const {
            assert(i1 <= i2); assert(i2 <= i3);
            auto i1_ex = m_gct_bp_support.excess(i1);
            auto i2_ex = m_gct_bp_support.excess(i2);
            auto i3_ex = m_gct_bp_support.excess(i3);
            //std::cout << i1_ex << " " << i2_ex << " " << i3_ex << std::endl;
            if(i1_ex < i2_ex) {
                if(i1_ex < i3_ex) return i1;
                else return i3;
            } else if(i1_ex > i2_ex) {
                if(i2_ex < i3_ex) return i2;
                else return i3;
            } else {
                if(i2_ex < i3_ex) return i2;
                else return i3;
            }
        }
        
        
    public:
        typedef typename bit_vector::size_type size_type;
        typedef typename bit_vector::size_type value_type;
        typedef t_bp_support                   bp_support_type;

        const bit_vector&      gct_bp         = m_gct_bp;
        const bp_support_type& gct_bp_support = m_gct_bp_support;

        //! Default constructor
        rmq_succinct_bp_fast() {}

        //! Constructor
        /*! \tparam t_rac A random access container.
         *  \param  v     Pointer to container object.
         */
        template<class t_rac>
        rmq_succinct_bp_fast(const t_rac* v=nullptr) {
            if (v != nullptr) {
                m_gct_bp = bit_vector(2*v->size()+2,0);
                construct_generalized_cartesian_tree_leftmost(v);
                m_gct_bp_support = bp_support_type(&m_gct_bp); 
                construct_sparse_table_over_sampled_excess_value_of_bp();
            }
        }

        //! Copy constructor
        rmq_succinct_bp_fast(const rmq_succinct_bp_fast& rm) {
            *this = rm;
        }

        //! Move constructor
        rmq_succinct_bp_fast(rmq_succinct_bp_fast&& rm) {
            *this = std::move(rm);
        }

        rmq_succinct_bp_fast& operator=(const rmq_succinct_bp_fast& rm) {
            if (this != &rm) {
                m_gct_bp = rm.m_gct_bp;
                m_gct_bp_support = rm.m_gct_bp_support;
                m_gct_bp_support.set_vector(&m_gct_bp);
                m_sparse_rmq = rm.m_sparse_rmq;
            }
            return *this;
        }

        rmq_succinct_bp_fast& operator=(rmq_succinct_bp_fast&& rm) {
            if (this != &rm) {
                m_gct_bp = std::move(rm.m_gct_bp);
                m_gct_bp_support = std::move(rm.m_gct_bp_support);
                m_gct_bp_support.set_vector(&m_gct_bp);   
                m_sparse_rmq = std::move(rm.m_sparse_rmq);
            }
            return *this;
        }

        void swap(rmq_succinct_bp_fast& rm) {
            m_gct_bp.swap(rm.m_gct_bp);
            util::swap_support(m_gct_bp_support, rm.m_gct_bp_support,
                               &m_gct_bp, &(rm.m_gct_bp));   
            m_sparse_rmq.swap(rm.m_sparse_rmq);
        }
        
        void print_bp() {
            size_t N = m_gct_bp.size();
            for(size_t i = 0; i < N; ++i) {
                std::cout << (m_gct_bp[i] ? '(' : ')');
            }
            std::cout << std::endl;
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
        size_type operator()(const size_type l, const size_type r)const {
            assert(l <= r); assert(r < size());
            size_type i     = m_gct_bp_support.select(l+2)-1;
            size_type j     = m_gct_bp_support.select(r+2);
            size_type sparse_i = i/t_block_size + (i % t_block_size != 0);
            size_type sparse_j = j/t_block_size-(j >= t_block_size);
            if(sparse_i >= sparse_j) {
                size_type rmq_e = m_gct_bp_support.rmq(i,j);
                return m_gct_bp_support.rank(rmq_e)-1;
            }
            size_type rmq_sparse = m_min_excess_idx[m_sparse_rmq(sparse_i,sparse_j)];
            size_type rmq_e1 = m_gct_bp_support.rmq(i,t_block_size*sparse_i-(sparse_i != 0 && i % t_block_size != 0));
            size_type rmq_e2 = m_gct_bp_support.rmq(t_block_size*(sparse_j+1),j);
            size_type rmq_e = min_excess(rmq_e1,rmq_sparse,rmq_e2);
            return m_gct_bp_support.rank(rmq_e)-1;
        }

        size_type size()const {
            return (m_gct_bp.size()-2)/2;
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_gct_bp.serialize(out, child, "gct_bp");
            written_bytes += m_gct_bp_support.serialize(out, child, "gct_bp_support");    
            written_bytes += m_sparse_rmq.serialize(out, child, "sparse_rmq"); 
            written_bytes += m_min_excess.serialize(out, child, "min_excess");             
            written_bytes += m_min_excess_idx.serialize(out, child, "min_excess_idx"); 
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            m_gct_bp.load(in);
            m_gct_bp_support.load(in, &m_gct_bp);
            //TODO: Load sparse table rmq
        }
};

} // end namespace sdsl
#endif
