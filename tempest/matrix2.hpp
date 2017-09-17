#pragma once
#ifndef TM_MATRIX2_HPP
#define TM_MATRIX2_HPP

/**
 *  @file matrix2.hpp
 *  @brief Defines tempest::matrix<T,2,2>.
 *  @author Toru Matsuoka
 *
 *  $Id: matrix2.hpp,v 1.7 2005/05/22 10:40:00      Exp $
 */

//-----------------------------------------------

//#include<cstddef>//std::size_t
//#include<cmath>
//#include<algorithm>
//#include<functional>
#include <cassert>

#include "matrix.hpp"
//#include"matrix_base.hpp"

namespace tempest
{
    
/**
 *  @class matrix<T,2,2>
 *  @brief matrix class template (Specific 2-dimentional)
 *
 *  @code
 *  tempest::matrix<double,2> m;    //2-dimentional matrix m
 *  @endcode
 *
 */

    template <class T>
    class matrix<T, 2, 2> : public matrix_base<matrix<T, 2, 2>, T, 2, 2>
    {
    public:
        //-----------------------------------------------
        //type defines
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef matrix<T, 2, 2> this_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        typedef T* iterator;
        typedef const T* const_iterator;

    public:
        static const size_type row_size = 2;
        static const size_type col_size = 2;
        static const size_type c_size = 2 * 2;

    public:
        //-----------------------------------------------
        //functions for iterator
        iterator begin() { return &(element[0][0]); }
        iterator end() { return &(element[0][0]) + c_size; }

        const_iterator begin() const { return &(element[0][0]); }
        const_iterator end() const { return &(element[0][0]) + c_size; }

        //-----------------------------------------------
        //constructors and destructor
        matrix() {}

        explicit matrix(
            T _m00, T _m01,
            T _m10, T _m11) : m00(_m00), m01(_m01),
                              m10(_m10), m11(_m11) {}

        matrix(const this_type& rhs) : m00(rhs.m00), m01(rhs.m01),
                                       m10(rhs.m10), m11(rhs.m11) {}

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
#endif

        template <class X>
        explicit matrix(const matrix<X, 2, 2>& rhs)
        {
            std::copy(rhs.begin(), rhs.end(), begin());
        }

        template <class X>
        explicit matrix(const X rhs[c_size])
        {
            std::copy(rhs, rhs + c_size, begin());
        }

        template <class X>
        explicit matrix(const X(&rhs)[row_size][col_size])
        {
            std::copy(&(rhs[0][0]), &(rhs[0][0]) + c_size, begin());
        }

#ifdef _MSC_VER
#pragma warning(pop)
#endif

        template <class Self, class Type>
        explicit matrix(const matrix_base<Self, Type, row_size, col_size>& rhs) : m00(static_cast<const Self&>(rhs[0][0])), m01(static_cast<const Self&>(rhs[0][1])),
                                                                                  m10(static_cast<const Self&>(rhs[1][0])), m11(static_cast<const Self&>(rhs[1][1]))
        {
        }

        ~matrix() {}

        //swap :-)
        //void swap(this_type& rhs){std::swap(m,rhs.m);}

        //-----------------------------------------------
        //inserters
        template <class X>
        this_type& operator=(const matrix<X, 2, 2>& rhs)
        {
            m00 = static_cast<T>(rhs[0][0]);
            m01 = static_cast<T>(rhs[0][1]);
            m10 = static_cast<T>(rhs[1][0]);
            m11 = static_cast<T>(rhs[1][1]);

            return *this;
        }

        template <class Self, class Type>
        this_type& operator=(const matrix_base<Self, Type, row_size, col_size>& rhs)
        {
            m00 = static_cast<T>(static_cast<const Self&>(rhs)[0][0]);
            m01 = static_cast<T>(static_cast<const Self&>(rhs)[0][1]);
            m10 = static_cast<T>(static_cast<const Self&>(rhs)[1][0]);
            m11 = static_cast<T>(static_cast<const Self&>(rhs)[1][1]);

            return *this;
        }

        //-----------------------------------------------
        //capacity
        size_type size() const { return c_size; }
        size_type max_size() const { return c_size; }
        bool empty() const { return false; }

        //-----------------------------------------------
        //operators
        this_type& negate()
        {
            m00 = -m00;
            m01 = -m01;
            m10 = -m10;
            m11 = -m11;
            return *this;
        }

//std::transform(begin(),end(),rhs.begin(),begin(),std::plus<T>());//
#define DECLARE_OP_EQUAL(OP)                     \
    this_type& operator OP(const this_type& rhs) \
    {                                            \
        m00 OP rhs.m00;                          \
        m01 OP rhs.m01;                          \
        m10 OP rhs.m10;                          \
        m11 OP rhs.m11;                          \
        return *this;                            \
    }

        DECLARE_OP_EQUAL(+= )
        DECLARE_OP_EQUAL(-= )

#undef DECLARE_OP_EQUAL

        //This functon is NOT recommended.
        this_type& operator*=(const this_type& rhs)
        {
            const T tp0 = m00;
            const T tp1 = m10;

            m00 = tp0 * rhs.m00 + m01 * rhs.m10;
            m01 = tp0 * rhs.m01 + m01 * rhs.m11;
            m10 = tp1 * rhs.m00 + m11 * rhs.m10;
            m11 = tp1 * rhs.m01 + m11 * rhs.m11;

            return *this;
        }

        //template<class X>
        this_type& operator*=(value_type rhs)
        {
            m00 *= rhs;
            m01 *= rhs;
            m10 *= rhs;
            m11 *= rhs;

            return *this;
        }

        //template<class X>
        this_type& operator/=(value_type rhs)
        {
            m00 /= rhs;
            m01 /= rhs;
            m10 /= rhs;
            m11 /= rhs;

            return *this;
        }

        T* operator[](size_type i)
        { //pointer what self is const.
            return element[i];
        }
        const T* operator[](size_type i) const
        { //pointer what self is const.
            return element[i];
        }

        //-----------------------------------------------
        // utilities
        this_type& transpose()
        {
            std::swap(m01, m10);
            return *this;
        }

        T det() const
        { //determinant
            return m00 * m11 - m01 * m10;
        }

        bool is_invertible() const
        {                              // nonsingular
            return (this->det() != 0); //|M| != 0
        }

        const char* debug() const { return "tempest::matrix<T,2,2>"; }

    private:
        union
        {
            struct
            {
                T m00, m01;
                T m10, m11;
            };
            T element[2][2];
        };
    };

    //-----------------------------------------------
    //Not a member!
    //-----------------------------------------------

    //-----------------------------------------------
    //operators

    template <class T> //inverter
    inline matrix<T, 2, 2> invert(const matrix<T, 2, 2>& rhs)
    {
        T det = rhs.det();

        assert(det != T());

        det = T(1) / det;

        return matrix<T, 2, 2>(
            rhs[1][1] * det, -rhs[0][1] * det,
            -rhs[1][0] * det, rhs[0][0] * det);
    }

    template <class T> //inverter
    inline matrix<T, 2, 2> operator~(const matrix<T, 2, 2>& rhs)
    {
        return invert(rhs);
    }

#ifdef minor
//#error ?
#else
    template <class T> //minor
    inline matrix<T, 2, 2> minor(const matrix<T, 2, 2>& rhs)
    {
        return matrix<T, 2, 2>(
            rhs[1][1], -rhs[0][1],
            -rhs[1][0], rhs[0][0]);
    }
#endif

    template <class T>
    matrix<T, 2, 2> operator*(const matrix<T, 2, 2>& lhs, const matrix<T, 2, 2>& rhs)
    {
        return matrix<T, 2, 2>(
            lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0],
            lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1],

            lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0],
            lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1]);
    }

    template <class T>
    void multiply(matrix<T, 2, 2>* out, const matrix<T, 2, 2>& lhs, const matrix<T, 2, 2>& rhs)
    {
        (*out)[0][0] = lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0];
        (*out)[0][1] = lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1];

        (*out)[1][0] = lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0];
        (*out)[1][1] = lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1];

        return;
    }

    //-----------------------------------------------
    // utilities
    template <class T>
    inline matrix<T, 2, 2> transpose(const matrix<T, 2, 2>& rhs)
    {
        return matrix<T, 2, 2>(
            rhs[0][0], rhs[1][0],
            rhs[0][1], rhs[1][1]);
    }

    template <class T>
    inline T det(const matrix<T, 2, 2>& rhs)
    {
        return rhs.det();
    }

    //--------------------------------------------------
    //compare

    template <class T>
    inline bool operator==(const matrix<T, 2, 2>& lhs, const matrix<T, 2, 2>& rhs)
    {
        return (lhs[0][0] == rhs[0][0]) && (lhs[0][1] == rhs[0][1]) &&
               (lhs[1][0] == rhs[1][0]) && (lhs[1][1] == rhs[1][1]);
    }

    template <class T>
    inline bool operator!=(const matrix<T, 2, 2>& lhs, const matrix<T, 2, 2>& rhs)
    {
        return !(lhs == rhs);
    }

} //end of namespace

//#include "matrix_functions.hpp"

#endif
