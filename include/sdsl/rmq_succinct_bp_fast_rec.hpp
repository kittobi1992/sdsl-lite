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
/*! \file rmq_succinct_bp_fast_rec.hpp
    \brief rmq_succinct_bp_fast_rec.hpp contains the class rmq_succinct_bp_fast which supports range minimum or range maximum queries on a random access container in constant time and \f$2 n+o(n) bits\f$ space.
    \author Tobias Heuer
*/
#ifndef INCLUDED_SDSL_RMQ_SUCCINCT_BP_FAST_REC
#define INCLUDED_SDSL_RMQ_SUCCINCT_BP_FAST_REC

#include <stack>
#include <limits>

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "bits.hpp"
#include "bp_support_sada.hpp"
#include "bp_support_algorithm.hpp"
#include "suffix_tree_helper.hpp"
#include "util.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint32_t t_block_size=512, bool fast = true, bool t_min = true,
         class t_rank = rank_support_v5<> >
class rmq_succinct_bp_fast_rec;

template<uint32_t t_block_size=512, bool fast = true, class t_rank = rank_support_v5<> >
struct range_maximum_bp_fast_rec {
    typedef rmq_succinct_bp_fast_rec<t_block_size, fast, false> type;
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
template<uint32_t t_block_size, bool fast, bool t_min, class t_rank>
class rmq_succinct_bp_fast_rec
{
        bit_vector                  m_gct_bp;         //!< A bit vector which contains the BP-GT of the input container.
        t_rank                      m_bp_rank; //!< Support structure for the BPS-GT
        rmq_succinct_bp_fast<128,fast,t_min,true>  m_rmq_recursive;
        int_vector<>                m_min_excess_idx;
        int_vector<>                m_min_excess; 
        int_vector<>                m_min_excess_idx64;
        bit_vector::value_type    m_max_depth;

        void copy(const rmq_succinct_bp_fast_rec& rm) {
            m_gct_bp = rm.m_gct_bp;
            m_bp_rank = rm.m_bp_rank;
            m_bp_rank.set_vector(&m_gct_bp);
            m_min_excess = rm.m_min_excess;
            m_min_excess_idx = rm.min_excess_idx;
            m_rmq_recursive = rm.m_rmq_recursive;
            m_min_excess_idx64 = rm.m_min_excess_idx64;
            m_max_depth = rm.m_max_depth;
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
            size_type bp_size = m_gct_bp.size();
            m_min_excess = int_vector<>(bp_size/t_block_size+1,0);
            m_min_excess_idx = int_vector<>(bp_size/t_block_size+1,0);
            for(size_t i = 0; i*t_block_size < bp_size; ++i) {
                uint64_t min_idx = fast_rmq_scan(i*t_block_size, std::min((i+1)*t_block_size - 1,bp_size-1));
                m_min_excess_idx[i] = min_idx-i*t_block_size;
                m_min_excess[i] = excess(min_idx);
            }
            m_rmq_recursive = rmq_succinct_bp_fast<128,fast,t_min,true>(&m_min_excess);
            util::bit_compress(m_min_excess);
            util::bit_compress(m_min_excess_idx);
        }
        
        void construct_minima_for_64_bit_blocks() {
            m_max_depth = std::numeric_limits<bit_vector::value_type>::max();
            size_type bp_size = m_gct_bp.size();
            m_min_excess_idx64 = int_vector<>(bp_size/64+1,0);
            size_type cur_excess = 0;
            for(size_t i = 0; i*64 < bp_size; ++i) {
                uint64_t min_idx = 0;
                size_type min_excess = std::numeric_limits<size_type>::max();
                for(size_t j = i*64; j <= std::min((i+1)*64 - 1,bp_size-1); ++j) {
                    if(m_gct_bp[j]) cur_excess++;
                    else cur_excess--;
                    if(cur_excess <= min_excess) {
                        min_idx = j;
                        min_excess = cur_excess;
                    }
                    if(cur_excess > m_max_depth) m_max_depth = cur_excess;
                }
                m_min_excess_idx64[i] = (min_idx-i*64)/8;
            }
            util::bit_compress(m_min_excess_idx64);
        }
        
        inline std::pair<bit_vector::size_type, int_vector<>::value_type> min_ex(const bit_vector::size_type i1, const bit_vector::size_type i2, const bit_vector::size_type i3,
                                            const int_vector<>::value_type i1_ex, const int_vector<>::value_type i2_ex, const int_vector<>::value_type i3_ex) const {
            assert(i1 <= i2); assert(i2 <= i3);
            //std::cout << i1_ex << " " << i2_ex << " " << i3_ex << std::endl;
            if(i1_ex > i2_ex) {
                if(i2_ex < i3_ex) return std::make_pair(i2,i2_ex);
                else return  std::make_pair(i3,i3_ex);;
            } else if(i1_ex < i2_ex) {
                if(i1_ex < i3_ex) return  std::make_pair(i1,i1_ex);
                else return std::make_pair(i3,i3_ex);
            } else {
                if(i2_ex < i3_ex) return  std::make_pair(i2,i2_ex);
                else return  std::make_pair(i3,i3_ex);;
            }
        }
        
        inline int_vector<>::value_type get_min_excess_idx(size_t i) const {
            return m_min_excess_idx[i] + i * t_block_size;
        }
        
        
        inline bit_vector::size_type get_min_pos_in_64bit_word(const uint64_t *bv_data, const bit_vector::size_type& block_id, const bit_vector::size_type& idx) const {
            typedef bit_vector::size_type size_type;
            size_type min_pos_byte_idx = idx+8*m_min_excess_idx64[block_id];
            size_type min_pos_idx = min_pos_byte_idx 
                                    + excess::data.min_pos_max[(((*(bv_data+(min_pos_byte_idx>>6)))>>(min_pos_byte_idx&0x3F))&0xFF)]; 
            return min_pos_idx;
        }
        
        inline bit_vector::difference_type get_min_excess_in_64bit_word(const uint64_t *bv_data, const bit_vector::size_type& block_id, const bit_vector::size_type& min_pos) const {
            value_type min_pos_data = *(bv_data+block_id) & bits::lo_set[min_pos];
            return (bits::cnt(min_pos_data)<<1) - min_pos;           
        }
        
        inline int_vector<>::value_type excess(size_t i) const {
            return (m_bp_rank(i+1)<<1)-(i+1);
        }
        
        inline bit_vector::size_type select(size_t i) const {
            typedef bit_vector::size_type size_type;
            size_type l = 0; size_type r = std::min(static_cast<size_type>(m_max_depth),static_cast<size_type>(2*(i-1))); //TODO: Calculate max depth of cartesian tree
            while(l < r) {
                size_type m = (l+r)/2 + ((l+r)&1);
                if(m_bp_rank(2*(i-1)-m+1) < i) r = m-1;
                else l = m;
            }
            return 2*(i-1)-l;
        }
        
        inline bit_vector::size_type naive_select(size_t j) const {
            typedef bit_vector::size_type size_type;
            size_type cnt = 0;
            for(size_t i = 0; i < m_gct_bp.size(); ++i) {
                if(m_gct_bp[i]) cnt++;
                if(cnt == j) return i;
            }
            return 0;
        }
        
        inline bit_vector::size_type fast_rmq_scan(const bit_vector::size_type l, const bit_vector::size_type r) const {
            typedef bit_vector::size_type size_type;
            typedef bit_vector::value_type value_type;
            typedef bit_vector::difference_type difference_type;
            
            
            const size_type l64 = ((l+64)/64)*64;
            const size_type r64 = (r/64)*64;
            difference_type cur_excess = 0;
            difference_type min_excess = 0;
            size_type min_excess_pos = l;
            const uint64_t* b = m_gct_bp.data();
            
//               std::cout << l << " " << l64 << " " << r64 << " " << r << " ("<<(l%64) << ","<< (l64%64) << "," << (r64%64)<<","<<(r%63) << ")" << std::endl;
            size_type min_pos_in_first_block = get_min_pos_in_64bit_word(b,l>>6,(l/64)*64);
            if(min_pos_in_first_block < l || min_pos_in_first_block > std::min(l64,r)) { 
                size_type l8 = (((l+1)+7)/8)*8;
                size_type r8 = (r/8)*8;
                
                //If l start inside an 8-bit value, we have determine the
                //excess value until the border of first 8-bit word inside (l,r)
                for(size_type i = l+1; i < std::min(l8,r+1); ++i) {
                    if(m_gct_bp[i]) ++cur_excess;
                    else {
                        --cur_excess;
                        if(cur_excess <= min_excess) {
                            min_excess_pos = i;
                            min_excess = cur_excess;
                        }
                    }
                }
                
                for(size_type i = l8; i < std::min(r8,l64); i += 8) {
                    uint8_t val = (((*(b+(i>>6)))>>(i&0x3F))&0xFF);
                    difference_type x = excess::data.min[val];
                    if(cur_excess+x <= min_excess) {
                        min_excess = cur_excess+x;
                        min_excess_pos = i + excess::data.min_pos_max[val];
                    }
                    cur_excess += excess::data.word_sum[val];
                }
                
                //Special case: If l starts in an 8-bit word and r end in the
                //following 8-bit word.
                for(size_type i = std::max(l8,r8); i < std::min(l64,r+1); ++i) {
                    if(m_gct_bp[i]) ++cur_excess;
                    else {
                        --cur_excess;
                        if(cur_excess <= min_excess) {
                            min_excess_pos = i;
                            min_excess = cur_excess;
                        }
                    }
                }
                
            }
            else {
                min_excess_pos = min_pos_in_first_block;
                min_excess = get_min_excess_in_64bit_word(b,l>>6,min_excess_pos-(l/64)*64+1);
                cur_excess += (bits::cnt(*(b+(l>>6)))<<1)-64;
            }
            
            
            //Scanning 64-bit words for min excess
            for(size_type i = l64; i < r64; i+=64) {
                size_type cur_block = i>>6;
                size_type min_pos_idx = get_min_pos_in_64bit_word(b,cur_block,i);
                difference_type tmp_min_excess = cur_excess + get_min_excess_in_64bit_word(b,cur_block,min_pos_idx-i+1);
                if(tmp_min_excess <= min_excess) {
                     min_excess_pos = min_pos_idx;
                     min_excess = tmp_min_excess;
                }
                cur_excess += (bits::cnt(*(b+cur_block))<<1)-64;
            }
            
            
            size_type last_block = std::max(l64,r64);
            if(last_block <= r) {
                size_type min_pos_in_last_block = get_min_pos_in_64bit_word(b,last_block>>6,(last_block/64)*64);
                if(min_pos_in_last_block > r) {
                
                    size_type r8 = (r/8)*8;
                
                    for(size_type i = last_block; i < r8; i += 8) {
                        uint8_t val = (((*(b+(i>>6)))>>(i&0x3F))&0xFF);
                        difference_type x = excess::data.min[val];
                        if(cur_excess+x <= min_excess) {
                            min_excess = cur_excess+x;
                            min_excess_pos = i + excess::data.min_pos_max[val];
                        }
                        cur_excess += excess::data.word_sum[val];
                    }
                
                    //If r ends inside an 8-bit value, we have determine the
                    //excess value until the border of last 8-bit word inside (l,r)
                    for(size_type i = std::max(r8,last_block); i < r+1; ++i) {
                        if(m_gct_bp[i]) ++cur_excess;
                        else {
                            --cur_excess;
                            if(cur_excess <= min_excess) {
                                min_excess_pos = i;
                                min_excess = cur_excess;
                            }
                        }
                    }
                
                }
                else {
                    difference_type tmp_min_excess = cur_excess + get_min_excess_in_64bit_word(b,last_block>>6,min_pos_in_last_block-(last_block/64)*64+1);
                    if(tmp_min_excess <= min_excess) {
                        min_excess_pos = min_pos_in_last_block;
                        min_excess = tmp_min_excess;
                    }
                }
            }
            
            
            /*size_type real_min_pos = m_gct_bp_support.rmq(l,r);
            if(real_min_pos != min_excess_pos) {
                std::cout << "Real Min Pos: " << real_min_pos << " vs. Current Min Pos: " << min_excess_pos << std::endl;
                abort();
            }*/
            
            
            return min_excess_pos;
        }
        
        
    public:
        typedef typename bit_vector::size_type size_type;
        typedef typename bit_vector::size_type value_type;
        typedef typename int_vector<>::value_type i_value_type;

        const bit_vector&                  gct_bp         = m_gct_bp;
        const t_rank&                      bp_rank        = m_bp_rank;
        rmq_succinct_bp_fast<128,fast,t_min,true>&  rmq_recursive = m_rmq_recursive;
        const int_vector<>&                min_excess     = m_min_excess;
        const int_vector<>&                min_excess_idx = m_min_excess_idx;
        const int_vector<>&                min_excess_idx64 = m_min_excess_idx64;

        //! Default constructor
        rmq_succinct_bp_fast_rec() {}

        //! Constructor
        /*! \tparam t_rac A random access container.
         *  \param  v     Pointer to container object.
         */
        template<class t_rac>
        rmq_succinct_bp_fast_rec(const t_rac* v=nullptr) {
            if (v != nullptr) {
                m_gct_bp = bit_vector(2*v->size()+2,0);
                construct_generalized_cartesian_tree_leftmost(v);
                construct_minima_for_64_bit_blocks();
                m_bp_rank = t_rank(&m_gct_bp); 
                construct_sparse_table_over_sampled_excess_value_of_bp();
            }
        }

        //! Copy constructor
        rmq_succinct_bp_fast_rec(const rmq_succinct_bp_fast_rec& rm) {
            *this = rm;
        }

        //! Move constructor
        rmq_succinct_bp_fast_rec(rmq_succinct_bp_fast_rec&& rm) {
            *this = std::move(rm);
        }

        rmq_succinct_bp_fast_rec& operator=(const rmq_succinct_bp_fast_rec& rm) {
            if (this != &rm) {
                m_gct_bp = rm.m_gct_bp;
                m_bp_rank = rm.m_bp_rank;
                m_bp_rank.set_vector(&m_gct_bp);
                m_min_excess = rm.m_min_excess;
                m_min_excess_idx = rm.min_excess_idx;
                m_rmq_recursive = rm.m_rmq_recursive;
                m_min_excess_idx64 = rm.min_excess_idx64;
                m_max_depth = rm.m_max_depth;
            }
            return *this;
        }

        rmq_succinct_bp_fast_rec& operator=(rmq_succinct_bp_fast_rec&& rm) {
            if (this != &rm) {
                m_gct_bp = std::move(rm.m_gct_bp);
                m_bp_rank = std::move(rm.m_bp_rank);
                m_bp_rank.set_vector(&m_gct_bp);   
                m_min_excess = std::move(rm.m_min_excess);
                m_min_excess_idx = std::move(rm.m_min_excess_idx);
                m_rmq_recursive = std::move(rm.m_rmq_recursive);
                m_min_excess_idx64 = std::move(rm.m_min_excess_idx64);
                m_max_depth = std::move(rm.m_max_depth);
            }
            return *this;
        }

        void swap(rmq_succinct_bp_fast_rec& rm) {
            m_gct_bp.swap(rm.m_gct_bp);
            util::swap_support(m_bp_rank, rm.m_bp_rank,
                               &m_gct_bp, &(rm.m_gct_bp));   
            m_min_excess.swap(rm.m_min_excess);
            m_min_excess_idx.swap(rm.m_min_excess_idx);
            m_rmq_recursive.swap(rm.m_rmq_recursive);
            m_min_excess_idx64.swap(rm.m_min_excess_idx64);
            m_max_depth = rm.m_max_depth;
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
        /*size_type operator()(const size_type l, const size_type r)const {
            assert(l <= r); assert(r < size());
            size_type i     = select(l+2)-1;
            size_type j     = select(r+2);
            size_type sparse_i = (i+t_block_size-1)/t_block_size;
            size_type sparse_j = j/t_block_size;
            if(sparse_i >= sparse_j) {
                sparse_i = sparse_i-(i != 0);
                size_type rmq_e = i;
                if(sparse_i == sparse_j) rmq_e = get_min_excess_idx(sparse_i);
                else {
                    size_type rmq_e1 = get_min_excess_idx(sparse_i);
                    size_type rmq_e2 = get_min_excess_idx(sparse_j);
                    rmq_e = (m_min_excess[sparse_i] < m_min_excess[sparse_j] ? rmq_e1 : rmq_e2);   
                }
                if(rmq_e < i || rmq_e > j) {
                    rmq_e = fast_rmq_scan(i,j);
                }
                return m_bp_rank(rmq_e+1)-1;
            }
            else {
                size_type rmq_sparse_idx = m_rmq_recursive(sparse_i,sparse_j-1);
                size_type rmq_sparse = get_min_excess_idx(rmq_sparse_idx);
                i_value_type rmq_sparse_ex = m_min_excess[rmq_sparse_idx];
                
                size_type rmq_e1 = get_min_excess_idx(sparse_i-(i != 0));
                i_value_type rmq_e1_ex = m_min_excess[sparse_i-(i != 0)];
                if(rmq_e1_ex < rmq_sparse_ex && rmq_e1 < i) {
                    rmq_e1 = fast_rmq_scan(i,t_block_size*sparse_i);
                    rmq_e1_ex = excess(rmq_e1);
                }
                
                size_type rmq_e2 = get_min_excess_idx(sparse_j);
                i_value_type rmq_e2_ex = m_min_excess[sparse_j];
                if(rmq_e2_ex <= rmq_sparse_ex && rmq_e2 > j) {
                    rmq_e2 = fast_rmq_scan(t_block_size*sparse_j,j);
                    rmq_e2_ex = excess(rmq_e2);
                }
                auto rmq_min = min_ex(rmq_e1,rmq_sparse,rmq_e2,
                                      rmq_e1_ex,rmq_sparse_ex,rmq_e2_ex);
                size_type rmq_e = rmq_min.first;
                int_vector<>::value_type rmq_e_ex = rmq_min.second;
                return ((rmq_e_ex+rmq_e)>>1);      
            }
        }*/
        
        std::pair<size_type,bool> operator()(const size_type l, const size_type r)const {
            assert(l <= r); assert(r < size());
            bool rmq_scan = false;
            size_type i     = select(l+2)-1;
            size_type j     = select(r+2);
            size_type sparse_i = (i+t_block_size-1)/t_block_size;
            size_type sparse_j = j/t_block_size;
            if(sparse_i >= sparse_j) {
                sparse_i = sparse_i-(i != 0);
                size_type rmq_e = i;
                if(sparse_i == sparse_j) rmq_e = get_min_excess_idx(sparse_i);
                else {
                    size_type rmq_e1 = get_min_excess_idx(sparse_i);
                    size_type rmq_e2 = get_min_excess_idx(sparse_j);
                    rmq_e = (m_min_excess[sparse_i] < m_min_excess[sparse_j] ? rmq_e1 : rmq_e2);   
                }
                if(rmq_e < i || rmq_e > j) {
                    rmq_scan = true;
                    rmq_e = fast_rmq_scan(i,j);
                }
                return std::make_pair(m_bp_rank(rmq_e+1)-1,rmq_scan);
            }
            else {
                size_type rmq_sparse_idx = m_rmq_recursive(sparse_i,sparse_j-1).first;
                size_type rmq_sparse = get_min_excess_idx(rmq_sparse_idx+1);
                i_value_type rmq_sparse_ex = m_min_excess[rmq_sparse_idx+1];
                
                size_type rmq_e1 = get_min_excess_idx(sparse_i-(i != 0));
                i_value_type rmq_e1_ex = m_min_excess[sparse_i-(i != 0)];
                if(rmq_e1_ex < rmq_sparse_ex && rmq_e1 < i) {
                    rmq_scan = true;
                    rmq_e1 = fast_rmq_scan(i,t_block_size*sparse_i);
                    rmq_e1_ex = excess(rmq_e1);
                }
                
                size_type rmq_e2 = get_min_excess_idx(sparse_j);
                i_value_type rmq_e2_ex = m_min_excess[sparse_j];
                if(rmq_e2_ex <= rmq_sparse_ex && rmq_e2 > j) {
                    rmq_scan = true;
                    rmq_e2 = fast_rmq_scan(t_block_size*sparse_j,j);
                    rmq_e2_ex = excess(rmq_e2);
                }
                
                /*size_t real = 0;
                size_t real_ex = std::numeric_limits<size_t>::max();
                for(size_t i = sparse_i; i < sparse_j; ++i) {
                    if(m_min_excess[i] <= real_ex) {
                        real_ex = m_min_excess[i];
                        real = get_min_excess_idx(i);
                    }
                }
                std::cout << (sparse_i-(i != 0)) << " " << sparse_i << " " << sparse_j-1 << " " << sparse_j << std::endl;
                std::cout << rmq_e1 << " " << rmq_sparse << " " << real << " " << rmq_e2 << std::endl; 
                std::cout << rmq_e1_ex << " " << rmq_sparse_ex << " " << real_ex << " " << rmq_e2_ex << std::endl; */
                
                auto rmq_min = min_ex(rmq_e1,rmq_sparse,rmq_e2,
                                      rmq_e1_ex,rmq_sparse_ex,rmq_e2_ex);
                size_type rmq_e = rmq_min.first;
                size_type rmq_e_ex = rmq_min.second;
                return std::make_pair((rmq_e_ex+rmq_e)/2-1,rmq_scan);      
            }
        }

        size_type size()const {
            return (m_gct_bp.size()-2)/2;
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_gct_bp.serialize(out, child, "gct_bp");
            written_bytes += m_bp_rank.serialize(out, child, "bp_rank");    
            written_bytes += m_min_excess.serialize(out, child, "min_excess");             
            written_bytes += m_min_excess_idx.serialize(out, child, "min_excess_idx"); 
            written_bytes += m_rmq_recursive.serialize(out, child, "rmq_recursive");
            written_bytes += m_min_excess_idx64.serialize(out, child, "min_excess_idx64"); 
            written_bytes += write_member(m_max_depth, out, child, "max_depth");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            m_gct_bp.load(in);
            m_bp_rank.load(in, &m_gct_bp);
            m_min_excess.load(in);
            m_min_excess_idx.load(in);
            m_rmq_recursive.load(in);
            m_min_excess_idx64.load(in);
            read_member(m_max_depth,in);
        }
};

} // end namespace sdsl
#endif
