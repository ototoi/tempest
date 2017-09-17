#pragma once
#ifndef TM_VECTOR2_HPP
#define TM_VECTOR2_HPP

/** 
 *  @file vector2.hpp
 *  @brief Defines tempest::vector<T,2>.
 *  @author Toru Matsuoka
 * 
 *  $Id: vector2.hpp,v 1.7 2005/03/17 13:00:40      Exp $
 */

#include <cstddef> //std::size_t
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "vector.hpp" //For static polimorphizm.

namespace tempest
{
    /** 
     *  @class vector<T,2>
     *  @brief vector class template (Specific 2-dimentional)
     *  
     *  @code
     *  tempest::vector<double,2> v;    
     *  @endcode
     *    
     */
    template <class T>
    class vector<T, 2> : public vector_base<vector<T, 2>, T, 2>
    {

    public:
        //-----------------------------------------------
        //type defines
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef vector<T, 2> this_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        typedef T* iterator;
        typedef const T* const_iterator;

    public:
        static const size_type c_size = 2; //container size

    public:
        //-----------------------------------------------
        //functions for iterator
        iterator begin() { return element; }
        iterator end() { return element + c_size; }
        const_iterator begin() const { return element; }
        const_iterator end() const { return element + c_size; }

        //-----------------------------------------------
        //front&back
        /**
     * @name front&back
     */
        //@{
        value_type front() { return element[0]; }
        value_type back() { return element[c_size - 1]; }
        value_type front() const { return element[0]; }
        value_type back() const { return element[c_size - 1]; }
        //@}

        //-----------------------------------------------
        //capacity
        size_type size() const { return c_size; }
        size_type max_size() const { return c_size; }
        bool empty() const { return false; }

#define SC_(a) static_cast<T>(a)
        //-----------------------------------------------
        //constructors and destructor
        vector() {}

        explicit vector(value_type _x) : m0(_x), m1(_x) {}

        vector(value_type _x, value_type _y) : m0(_x), m1(_y) {}

        vector(const this_type& rhs) : m0(rhs.m0), m1(rhs.m1) {}

        template <class X>
        explicit vector(const X rhs[c_size]) : m0(SC_(rhs[0])), m1(SC_(rhs[1]))
        {
        }

        template <class X>
        explicit vector(const vector<X, 2>& rhs) : m0(SC_(rhs[0])), m1(SC_(rhs[1]))
        {
        }

        template <class Self, class X>
        explicit vector(const vector_base<Self, X, 2>& rhs) : m0(SC_(static_cast<const Self&>(rhs)[0])),
                                                              m1(SC_(static_cast<const Self&>(rhs)[1]))
        {
        }

        ~vector() {}
        //
        //void swap(this_type& rhs){std::swap(*this,rhs);}

        //-----------------------------------------------
        //inserters
        this_type& operator=(const this_type& rhs)
        {
            m0 = rhs.m0;
            m1 = rhs.m1;
            return *this;
        }

        template <class X>
        this_type& operator=(const vector<X, 2>& rhs)
        {
            m0 = SC_(rhs[0]);
            m1 = SC_(rhs[1]);
            return *this;
        }

        template <class Self, class X>
        this_type& operator=(const vector_base<Self, X, 2>& rhs)
        {
            m0 = SC_(static_cast<const Self&>(rhs)[0]);
            m1 = SC_(static_cast<const Self&>(rhs)[1]);
            return *this;
        }

        template <class IT>
        void assign(IT start, IT end)
        {
            using namespace std;

            assert(distance(start, end) <= c_size); //debug
            copy(start, end, begin());
        }

        void assign(size_type num, value_type val)
        {
            size_t sz = (num < c_size) ? num : c_size;
            for (size_t i = 0; i < sz; i++)
                element[i] = val;
        }

#undef SC_

        //-----------------------------------------------
        //operators

        //--------------------------------
        this_type& negate()
        {
            m0 = -m0;
            m1 = -m1;
            return *this;
        }

#define DECLARE_OP_EQUAL(OP)                            \
    \
template<class X>                                       \
        this_type& operator OP(const vector<X, 2>& rhs) \
    {                                                   \
        m0 OP rhs.m0;                                   \
        m1 OP rhs.m1;                                   \
        return *this;                                   \
    }

        DECLARE_OP_EQUAL(+=)
        DECLARE_OP_EQUAL(-=)
        DECLARE_OP_EQUAL(*=)
        DECLARE_OP_EQUAL(/=)

#undef DECLARE_OP_EQUAL

        //--------------------------------
        template <class X>
        this_type& operator*=(const X rhs)
        {
            m0 *= rhs;
            m1 *= rhs;

            return *this;
        }
        template <class X>
        this_type& operator/=(const X rhs)
        {
            m0 /= rhs;
            m1 /= rhs;

            return *this;
        }
/*
    template<class X>
    this_type& operator*= (X rhs){
        m0 = static_cast<value_type>(m0 * rhs);
        m1 = static_cast<value_type>(m1 * rhs);
        return *this;
    }
    template<class X>
    this_type& operator/= (X rhs){
        m0 = static_cast<value_type>(m0 / rhs);
        m1 = static_cast<value_type>(m1 / rhs);
        return *this;
    }
*/

#undef DECLARE_OP_EQUAL_SCALAR

        //--------------------------------

        T& operator[](size_type i)
        {
            return element[i];
        }

        value_type operator[](size_type i) const
        {
            return element[i];
        }

        T& at(size_type i)
        {
            if (c_size <= i)
            {
                throw std::out_of_range(debug());
            }
            return element[i];
        }
        const T& at(size_type i) const
        {
            if (c_size <= i)
            {
                throw std::out_of_range(debug());
            }
            return element[i];
        }

        //-----------------------------------------------
        //utilities
        value_type length() const
        {
            using namespace std;
            return sqrt(sqr_length());
        }

        value_type sqr_length() const
        {
            return (
                (m0 * m0) + (m1 * m1));
        }

        value_type sum() const
        {
            return m0 + m1;
        }

        this_type& normalize()
        {
            using namespace std;

            value_type length = sqr_length(); //||V||^2

            length = T(1) / sqrt(length);
            m0 *= length;
            m1 *= length;

            return *this;
        }

        //-----------------------------------------------
        //debug code
    public:
        const char* debug() const { return "tempest::vector<T,2>"; }

    private:
        /** @union */
        union {
            /** @struct */
            struct
            {
                T m0, m1;
            };
            T element[2];
        };

    }; //End of class.

    //-----------------------------------------------
    //Not a member!
    //-----------------------------------------------

    //-----------------------------------------------
    //operators

    //See vector.hpp!

    //-----------------------------------------------
    //utility functions

    //normalize
    //length
    //sqr_length

    template <class T>
    inline T dot(const vector<T, 2>& lhs, const vector<T, 2>& rhs)
    {
        return (lhs[0] * rhs[0]) + (lhs[1] * rhs[1]);
    }

    //-----------------------------------------------
    //compare
    template <class T>
    inline bool operator==(const vector<T, 2>& lhs, const vector<T, 2>& rhs)
    {
        return (lhs[0] == rhs[0]) && (lhs[1] == rhs[1]);
    }

    //-----------------------------------------------
    //
} //End of namespace.

#include "vector_functions.hpp"

#endif
