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
/*! \file rank_select_support_bp.hpp
 *    \brief TODO
 *    \author Tobias Heuer
 */
#ifndef INCLUDED_SDSL_RANK_SELECT_SUPPORT_BP
#define INCLUDED_SDSL_RANK_SELECT_SUPPORT_BP

#include <stack>
#include <limits>

#include "rank_support_v5.hpp"
#include "int_vector.hpp"
#include "bits.hpp"
#include "util.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{
    
    
    //! TODO
    /*!
     * 
     * 
     */
    template<uint32_t t_sample_size = 1024>
    class rank_select_support_bp
    {
        
        const bit_vector*           m_v;
        rank_support_v5<>           m_bp_rank;
        int_vector<>                m_select_sample;
        bit_vector::value_type      m_max_excess;
        
        
        void copy(const rank_select_support_bp& rm) {
            m_v = rm.m_v;
            m_bp_rank = rm.m_bp_rank;
            m_bp_rank.set_vector(m_v);
            m_select_sample = rm.m_select_sample;
            m_max_excess = rm.m_max_excess;
        }
        
        
    private:
        
        void calculate_maximum_excess_value() {
            m_max_excess = std::numeric_limits<bit_vector::value_type>::min();
            bit_vector::value_type cur_excess = 0;
            for(size_t i = 0; i < m_v->size(); ++i) {
                if((*m_v)[i]) cur_excess++;
                else cur_excess--;
                if(cur_excess > m_max_excess) m_max_excess = cur_excess;
            }
            if(m_v->size() <= 2) m_max_excess = 1;
        }        
        
        void generate_select_sample() {
            size_t N = m_v->size()/2;
            m_select_sample = int_vector<>(N/t_sample_size+2,0);
            for(size_t i = 1; i < N; i += t_sample_size) {
                m_select_sample[i/t_sample_size] = select(i,false);
            }
            m_select_sample[N/t_sample_size+1] = select(N,false);
            util::bit_compress(m_select_sample);
        }
                                                                                 
                                                                                 
    public:
        typedef typename bit_vector::size_type size_type;
        typedef typename bit_vector::value_type value_type;
        typedef rank_support_trait<1,1>  trait_type;
        
        //! Constructor
        /*
         *  \param  v     Pointer to a bit vector
         */
        rank_select_support_bp(const bit_vector* v=nullptr) {
            if (v != nullptr) {
                set_vector(v);
                m_bp_rank = rank_support_v5<>(v);
                calculate_maximum_excess_value();
//                 std::cout << m_max_excess << std::endl;
                if(m_max_excess > t_sample_size) {
                    generate_select_sample();
                }
            }
        }
        
        //! Copy constructor
        rank_select_support_bp(const rank_select_support_bp& rm) {
            *this = rm;
        }
        
        //! Move constructor
        rank_select_support_bp(rank_select_support_bp&& rm) {
            *this = std::move(rm);
        }
        
        rank_select_support_bp& operator=(const rank_select_support_bp& rm) {
            if (this != &rm) {
                m_v = rm.m_v;
                m_bp_rank = rm.m_bp_rank;
                m_bp_rank.set_vector(m_v);
                m_select_sample = rm.m_select_sample;
                m_max_excess = rm.m_max_excess;
            }
            return *this;
        }
        
        rank_select_support_bp& operator=(rank_select_support_bp&& rm) {
            if (this != &rm) {
                m_v = std::move(rm.m_v);
                m_bp_rank = std::move(rm.m_bp_rank);
                m_bp_rank.set_vector(m_v);
                m_select_sample = std::move(rm.m_select_sample);
                m_max_excess = std::move(rm.m_max_excess);
            }
            return *this;
        }
        
        inline size_type excess(size_type idx) const {
            return (m_bp_rank(idx+1)<<1)-(idx+1);
        }
        
        inline size_type rank(size_type idx) const {
            return m_bp_rank.rank(idx);
        }
        
        /*inline size_type select(size_type idx, bool use_select_samples=true) const {
            std::cout << "Result: select("<<idx<<") = " << select2(idx) << std::endl;
            size_type l = 2*(idx-1)-std::min(static_cast<size_type>(m_max_excess),static_cast<size_type>(2*(idx-1))); size_type r = 2*(idx-1); 
            
            std::cout << "Start [l,r] = [" << l << "," << r << "]" << std::endl;
            
            if(use_select_samples && m_max_excess > t_sample_size) {
                size_t sample_idx = (idx-1)/t_sample_size;
                l = std::max(l,m_select_sample[sample_idx]);
                r = std::min(r,m_select_sample[sample_idx+1]);
            }
            std::cout << "Sample adaption [l,r] = [" << l << "," << r << "]" << std::endl;
            const uint64_t* p = m_bp_rank.m_basic_block.data()
                                + ((l>>10)&0xFFFFFFFFFFFFFFFEULL);// (idx/2048)*2
            //                     ( prefix sum of the 6x64bit blocks | (idx%2048)/(64*6) )
            size_type result = *p
                                + ((*(p+1)>>(60-12*((l&0x7FF)/(64*6))))&0x7FFULL)
                                + trait_type::word_rank(m_v->data(), l);
                                
            std::cout << result << " " << m_bp_rank(l) << std::endl;
            
            std::cout << "----------------------------------------" << std::endl;
            
            return 0;
        }*/
        
        inline size_type select(size_type idx, bool use_select_samples=true) const {
            size_type l = 2*(idx-1)-std::min(static_cast<size_type>(m_max_excess),static_cast<size_type>(2*(idx-1))); size_type r = 2*(idx-1); 
            
            if(use_select_samples && m_max_excess > t_sample_size) {
                size_t sample_idx = (idx-1)/t_sample_size;
                l = std::max(l,m_select_sample[sample_idx]);
                r = std::min(r,m_select_sample[sample_idx+1]);
            }
            
            while(r - l >= 64) {
                size_type m = (l+r)/2;
                if(m_bp_rank.rank(m) < idx) l = m+1;
                else r = m;
            }
            
            size_type s_pos = l;
            size_type s_rank = m_bp_rank(s_pos);
            value_type data = m_v->get_int(s_pos,r-l);
            r = r - l + 1; l = 0;
            while(l < r) {
                size_type m = (l+r)/2;
                if(s_rank + bits::cnt(data & bits::lo_set[m]) < idx) l = m+1;
                else r = m;
            }
            return s_pos + l - 1;
        }
        
        
        void swap(rank_select_support_bp& rm) {
            m_bp_rank.swap(rm.m_bp_rank); 
            m_select_sample.swap(rm.m_select_sample);
            m_max_excess = rm.m_max_excess;
        }
        
        
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_bp_rank.serialize(out, child, "bp_rank");
            written_bytes += m_select_sample.serialize(out, child, "select_sample");
            written_bytes += write_member(m_max_excess, out, child, "max_depth");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
        
        void load(std::istream& in, const bit_vector* v=nullptr) {
            m_bp_rank.load(in,v);
            m_select_sample.load(in);
            read_member(m_max_excess,in);
            set_vector(v);
        }
        
        void set_vector(const bit_vector* v=nullptr) {
            m_v = v;
            m_bp_rank.set_vector(m_v);
        }
        
    };
    
} // end namespace sdsl
#endif
