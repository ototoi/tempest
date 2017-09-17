#pragma once
#ifndef TM_VECTOR_HPP
#define TM_VECTOR_HPP

/** 
 * @file vector.hpp
 * @brief Defines tempest::vector<T,Sz>.
 *
 * $RCSfile: vector.hpp,v $
 * $Date: 2005/03/17 13:00:40 $
 * $Author: Toru Matsuoka $
 */

//-----------------------------------------------
#include <cstddef>   //std::size_t
#include <cmath>     //std::sqrt
#include <algorithm> //copy
#include <stdexcept>
#include <cassert>

#include "vector_base.hpp" //For static polimorphizm.

//_TEMPEST_USE_META_
//_TEMPEST_USE_FASTMATH_

namespace tempest
{

    /** 
     *  @brief Vector class template. 
     * 
     *  @code
     *  tempest::vector<double,3> v;    //3D vector class  v
     *  @endcode
     *
     */
    template <class T, std::size_t Sz>
    class vector : public vector_base<vector<T, Sz>, T, Sz>
    {
        //template<class,std::size_t> friend class vector
    public:
        //-----------------------------------------------
        //type define

        typedef T value_type;             /**< type of element*/
        typedef T& reference;             /**< reference*/
        typedef const T& const_reference; /**< const reference*/
        typedef vector<T, Sz> this_type;  /**< self*/

        typedef std::size_t size_type;          /**< type of element size*/
        typedef std::ptrdiff_t difference_type; /**< difference type of pointer*/

        typedef T* iterator;             /**< iterator*/
        typedef const T* const_iterator; /**< const iterator*/

    public:
        static const size_type c_size = Sz; //container size

    public:
        //-----------------------------------------------
        //function for iterator
        /** 
         * @name function for iterator
         */
        //@{
        iterator begin() { return element.begin(); }             /**< Return iterator for first element. */
        iterator end() { return element.end(); }                 /**< Return iterator for last element.  */
        const_iterator begin() const { return element.begin(); } /**< Return iterator for first element. */
        const_iterator end() const { return element.end(); }     /**< Return iterator for last element.  */
        //@}

        //-----------------------------------------------
        //front&back
        /**
         * @name front&back
         */
        //@{
        value_type front() { return element[0]; }
        value_type back() { return element[Sz - 1]; }
        const_reference front() const { return element[0]; }
        const_reference back() const { return element[Sz - 1]; }
        //@}

        //-----------------------------------------------
        //capacity
        /** 
         * @name size & capacity
         */
        //@{
        size_type size() const { return Sz; }     /**< Return size.           */
        size_type max_size() const { return Sz; } /**< Return size.           */
        bool empty() const { return false; }      /**< Return true always.    */
        //@}

        //-----------------------------------------------
        //constructor & destructor
        /** 
         * @name constructor & destructor
         */
        //@{
        vector() {}

        vector(const vector<T, Sz>& rhs)
            : element(rhs.element) //std::copy(rhs.begin(),rhs.end(),begin());
        {
            //detail::vector_hpp::node_loop<Sz-1>::eval(element,rhs);
            //node_loop<Sz-1>::eval(element,rhs);
        }

        template <class X>
        explicit vector(const vector<X, Sz>& rhs)
        {
            std::copy(rhs.begin(), rhs.end(), begin());
        }

        template <class X>
        explicit vector(const X* rhs)
        {
            std::copy(&(rhs[0]), &(rhs[Sz]), begin());
        }

        template <class Self, class X>
        explicit vector(const vector_base<Self, X, Sz>& rhs)
        {
            std::copy(static_cast<const Self&>(rhs).begin(), static_cast<const Self&>(rhs).end(), begin());
        }

        ~vector() {}
        //@}

        //swap :-)
        //void swap(this_type& rhs){std::swap(m,rhs.m);}

        //-----------------------------------------------
        //substitutioner
        /** 
         * @name substitutioner
         */
        //@{
        this_type& operator=(const this_type& rhs)
        {
            element = rhs.element;
            return *this;
        }
        template <class X>
        this_type& operator=(const vector<X, Sz>& rhs)
        {
            std::copy(rhs.begin(), rhs.end(), begin());
            return *this;
        }

        template <class Self, class X>
        this_type& operator=(const vector_base<Self, X, Sz>& rhs)
        {
            std::copy(static_cast<const Self&>(rhs).begin(), static_cast<const Self&>(rhs).end(), begin());
            return *this;
        }

        template <class IT>
        void assign(IT start, IT end)
        {
            assert(std::distance(start, end) <= c_size); //debug

            std::copy(start, end, begin());
        }

        void assign(size_type num, T val)
        {
            std::fill_n(begin(), (num < c_size) ? num : c_size, val);
        }

        //@}

        //-----------------------------------------------
        //operator
        /** 
         * @name operator
         */
        //@{

        //--------------------------------
        this_type& negate()
        {
            for (size_type i = 0; i < Sz; i++)
            {
                element[i] = -element[i];
            }
            return *this;
        }

//--------------------------------
#define DECLARE_OPERATOR_EQ(OP)                     \
template<class X>                                   \
this_type& operator OP(const vector<X, Sz>& rhs)    \
{                                                   \
    for (size_type i = 0; i < Sz; i++)              \
    {                                               \
        element[i] OP rhs[i];                       \
    }                                               \
    return *this;                                   \
}

        DECLARE_OPERATOR_EQ(+=)
        DECLARE_OPERATOR_EQ(-=)
        DECLARE_OPERATOR_EQ(*=)
        DECLARE_OPERATOR_EQ(/=)

#undef DECLARE_OPERATOR_EQ

        //--------------------------------
        template <class X>
        this_type& operator*=(const X rhs)
        {
            for (size_type i = 0; i < Sz; i++)
            {
                element[i] *= rhs;
            }
            return *this;
        }

        template <class X>
        this_type& operator/=(const X rhs)
        {
            for (size_type i = 0; i < Sz; i++)
            {
                element[i] /= rhs;
            }
            return *this;
        }

        //--------------------------------
        //@}

        //-----------------------------------------------
        //indexer
        /** 
         * @name indexer
         */
        //@{
        reference operator[](size_type i)
        {
            return element[i];
        }

        const value_type operator[](size_type i) const
        {
            return element[i];
        }

        reference at(size_type i)
        {
            if (Sz <= i)
            {
                throw std::out_of_range(debug());
            }
            return element[i];
        }
        const_reference at(size_type i) const
        {
            if (Sz <= i)
            {
                throw std::out_of_range(debug());
            }
            return element[i];
        }
        //@}

        //-----------------------------------------------
        //utility
        /** 
         * @name utility
         */
        //@{
        value_type length() const
        {
            using namespace std;
            return sqrt(sqr_length());
        }

        value_type sqr_length() const
        {
            T temp = T();
            for (size_type i = 0; i < Sz; i++)
            {
                temp += element[i] * element[i];
            }
            return temp;
        }

        value_type sum() const
        {
            T temp = T();
            for (size_type i = 0; i < Sz; i++)
            {
                temp += element[i];
            }
            return temp;
        }

        this_type& normalize()
        {
            T length = sqr_length(); //||V||^2
            //if (length == T()) return *this;
            //if (length){

            length = T(1) / std::sqrt(length);

            for (size_type i = 0; i < Sz; i++)
            {
                element[i] *= length;
            }

            return *this;
            //}
        }
        //@}

        /** 
         * @name others
         */
        //@{
        const char* debug() const { return "tempest::vector<T,Sz>"; }
        //@}

    private:
        struct ___m
        {
            T m[c_size];

            T* begin() { return m; }
            T* end() { return m + c_size; }
            const T* begin() const { return m; }
            const T* end() const { return m + c_size; }
            T& operator[](size_type i) { return m[i]; }
            const T operator[](size_type i) const { return m[i]; }
        } element;

    }; //End of class.

} //End of namespace.

#include "vector_functions.hpp"

#endif
