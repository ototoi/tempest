#pragma once
#ifndef TM_MATRIX4_HPP
#define TM_MATRIX4_HPP

/**
 * @file matrix4.hpp
 * @brief Defines tempest::matrix<T,4,4>.
 *
 * $RCSfile: matrix4.hpp,v $
 * $Date: 2005/05/22 10:39:47 $
 * $Author: Toru Matsuoka $
 */

//-----------------------------------------------

//#include<cstddef>//std::size_t
//#include<cmath>
//#include<algorithm>
#include <cstring> //memcpy
#include <functional>
#include <cassert>

#include "matrix.hpp"
//#include"matrix_base.hpp"

namespace tempest
{

    template <class T>
    inline T det33(
        T m00, T m01, T m02,
        T m10, T m11, T m12,
        T m20, T m21, T m22)
    {
        return (
            m00 * m11 * m22 - m11 * m20 * m02 +
            m01 * m12 * m20 - m10 * m22 * m01 +
            m02 * m10 * m21 - m12 * m21 * m00);
    }

/**
 *  @class matrix<T,4,4>
 *  @brief matrix class template (Specific 4-dimentional)
 *  
 *  @code
 *  tempest::matrix<double,4> m;    //4-dimentional matrix m
 *  @endcode
 *
 */

    template <class T>
    class matrix<T, 4, 4> : public matrix_base<matrix<T, 4, 4>, T, 4, 4>
    {

    public:
        //-----------------------------------------------
        //type defines
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef matrix<T, 4, 4> this_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        typedef T* iterator;
        typedef const T* const_iterator;

    public:
        //-----------------------------------------------
        //size
        static const size_type row_size = 4;
        static const size_type col_size = 4;
        static const size_type c_size = 4 * 4; //16

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
            T _m00, T _m01, T _m02, T _m03,
            T _m10, T _m11, T _m12, T _m13,
            T _m20, T _m21, T _m22, T _m23,
            T _m30, T _m31, T _m32, T _m33) : m00(_m00), m01(_m01), m02(_m02), m03(_m03),
                                              m10(_m10), m11(_m11), m12(_m12), m13(_m13),
                                              m20(_m20), m21(_m21), m22(_m22), m23(_m23),
                                              m30(_m30), m31(_m31), m32(_m32), m33(_m33) {}

        matrix(const this_type& rhs)
            : m00(rhs.m00), m01(rhs.m01), m02(rhs.m02), m03(rhs.m03),
              m10(rhs.m10), m11(rhs.m11), m12(rhs.m12), m13(rhs.m13),
              m20(rhs.m20), m21(rhs.m21), m22(rhs.m22), m23(rhs.m23),
              m30(rhs.m30), m31(rhs.m31), m32(rhs.m32), m33(rhs.m33)
        {
        }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
#endif

        template <class X>
        explicit matrix(const matrix<X, 4, 4>& rhs)
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
            std::copy(&rhs[0][0], &rhs[0][0] + c_size, begin());
        }

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifndef __clang__
        template <class Self, class Type>
        explicit matrix(const matrix_base<Self, Type, row_size, col_size>& rhs) : m00(static_cast<const Self&>(rhs)[0][0]), m01(static_cast<const Self&>(rhs)[0][1]), m02(static_cast<const Self&>(rhs)[0][2]), m03(static_cast<const Self&>(rhs)[0][3]),
                                                                                  m10(static_cast<const Self&>(rhs)[1][0]), m11(static_cast<const Self&>(rhs)[1][0]), m12(static_cast<const Self&>(rhs)[1][2]), m13(static_cast<const Self&>(rhs)[1][3]),
                                                                                  m20(static_cast<const Self&>(rhs)[2][0]), m21(static_cast<const Self&>(rhs)[2][1]), m22(static_cast<const Self&>(rhs)[2][2]), m23(static_cast<const Self&>(rhs)[2][3]),
                                                                                  m30(static_cast<const Self&>(rhs)[3][0]), m21(static_cast<const Self&>(rhs)[3][1]), m32(static_cast<const Self&>(rhs)[3][2]), m33(static_cast<const Self&>(rhs)[3][3])
        {
        }
#else
        template <class Self, class Type>
        explicit matrix(const matrix_base<Self, Type, row_size, col_size>& rhs)
        {
            m00 = static_cast<const Self&>(rhs)[0][0];
            m01 = static_cast<const Self&>(rhs)[0][1];
            m02 = static_cast<const Self&>(rhs)[0][2];
            m03 = static_cast<const Self&>(rhs)[0][3];
            m10 = static_cast<const Self&>(rhs)[1][0];
            m11 = static_cast<const Self&>(rhs)[1][0];
            m12 = static_cast<const Self&>(rhs)[1][2];
            m13 = static_cast<const Self&>(rhs)[1][3];
            m20 = static_cast<const Self&>(rhs)[2][0];
            m21 = static_cast<const Self&>(rhs)[2][1];
            m22 = static_cast<const Self&>(rhs)[2][2];
            m23 = static_cast<const Self&>(rhs)[2][3];
            m30 = static_cast<const Self&>(rhs)[3][0];
            m21 = static_cast<const Self&>(rhs)[3][1];
            m32 = static_cast<const Self&>(rhs)[3][2];
            m33 = static_cast<const Self&>(rhs)[3][3];
        }
#endif

        ~matrix()
        {
        }

        //swap :-)
        //void swap(this_type& rhs){std::swap(m,rhs.m);}

        //-----------------------------------------------
        //inserters
        template <class X>
        this_type& operator=(const matrix<X, 4, 4>& rhs)
        {
            std::copy(rhs.begin(), rhs.end(), begin());
            return *this;
        }

        template <class Self, class Type>
        this_type& operator=(const matrix_base<Self, Type, row_size, col_size>& rhs)
        {
            std::copy(static_cast<const Self&>(rhs).begin(), static_cast<const Self&>(rhs).end(), begin());
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
            for (size_type i = 0; i < row_size; ++i)
                for (size_type j = 0; j < col_size; ++j)
                    element[i][j] = -element[i][j];

            return *this;
        }

        //template<class X>
        this_type& operator+=(const this_type& rhs)
        {
#define DECLARE_OP(i, j) m##i##j += rhs.m##i##j;
            DECLARE_OP(0, 0)
            DECLARE_OP(0, 1)
            DECLARE_OP(0, 2)
            DECLARE_OP(0, 3)
            DECLARE_OP(1, 0)
            DECLARE_OP(1, 1)
            DECLARE_OP(1, 2)
            DECLARE_OP(1, 3)
            DECLARE_OP(2, 0)
            DECLARE_OP(2, 1)
            DECLARE_OP(2, 2)
            DECLARE_OP(2, 3)
            DECLARE_OP(3, 0)
            DECLARE_OP(3, 1)
            DECLARE_OP(3, 2)
            DECLARE_OP(3, 3)
#undef DECLARE_OP

            //std::transform(begin(),end(),rhs.begin(),begin(),std::plus<T>());//
            return *this;
        }

        //template<class X>
        this_type& operator-=(const this_type& rhs)
        {
#define DECLARE_OP(i, j) m##i##j -= rhs.m##i##j;
            DECLARE_OP(0, 0)
            DECLARE_OP(0, 1)
            DECLARE_OP(0, 2)
            DECLARE_OP(0, 3)
            DECLARE_OP(1, 0)
            DECLARE_OP(1, 1)
            DECLARE_OP(1, 2)
            DECLARE_OP(1, 3)
            DECLARE_OP(2, 0)
            DECLARE_OP(2, 1)
            DECLARE_OP(2, 2)
            DECLARE_OP(2, 3)
            DECLARE_OP(3, 0)
            DECLARE_OP(3, 1)
            DECLARE_OP(3, 2)
            DECLARE_OP(3, 3)
#undef DECLARE_OP
            //std::transform(begin(),end(),rhs.begin(),begin(),std::minus<T>());//
            return *this;
        }

        //This functon is NOT recommended.
        //template<class X>
        this_type& operator*=(const matrix<T, 4, 4>& rhs)
        {
            this_type temp(*this);

            //std::fill_n( (*this).begin(), c_size, T() );
            fill_zero_n((*this).begin(), c_size);

            for (size_type i = 0; i < 4; ++i)
            {
                for (size_type k = 0; k < 4; k++)
                {
                    for (size_type j = 0; j < 4; ++j)
                    {
                        element[i][j] += static_cast<T>(temp[i][k] * rhs[k][j]);
                    }
                }
            }
            return *this;
        }

        //template<class X>
        this_type& operator*=(value_type rhs)
        {
            for (size_type i = 0; i < row_size; ++i)
                for (size_type j = 0; j < col_size; ++j)
                    element[i][j] *= rhs;

            return *this;
        }
        //template<class X>
        this_type& operator/=(value_type rhs)
        {
            for (size_type i = 0; i < row_size; ++i)
                for (size_type j = 0; j < col_size; ++j)
                    element[i][j] /= rhs;

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
            std::swap(m02, m20);
            std::swap(m03, m30);
            std::swap(m12, m21);
            std::swap(m13, m31);
            std::swap(m23, m32);
            return *this;
        }

        T det() const
        { //derminant
            return (
                m00 * det33(m11, m12, m13, m21, m22, m23, m31, m32, m33) - m10 * det33(m01, m02, m03, m21, m22, m23, m31, m32, m33) + m20 * det33(m01, m02, m03, m11, m12, m13, m31, m32, m33) - m30 * det33(m01, m02, m03, m11, m12, m13, m21, m22, m23));
        }

        bool is_invertible() const
        {                                // nonsingular
            return (this->det() != T()); //|M| != 0
        }

        const char* debug() const { return "tempest::matrix<T,4,4>"; }

    private:
        union
        {
            struct
            {
                T m00, m01, m02, m03;
                T m10, m11, m12, m13;
                T m20, m21, m22, m23;
                T m30, m31, m32, m33;
            };
            T element[4][4];
        };
    };

    //-----------------------------------------------
    //Not a member!
    //-----------------------------------------------

    //-----------------------------------------------
    //operators

    template <class T>
    inline matrix<T, 4, 4> operator*(const matrix<T, 4, 4>& lhs, const matrix<T, 4, 4>& rhs)
    {
        return matrix<T, 4, 4>(
            lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0] + lhs[0][2] * rhs[2][0] + lhs[0][3] * rhs[3][0],
            lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1] + lhs[0][2] * rhs[2][1] + lhs[0][3] * rhs[3][1],
            lhs[0][0] * rhs[0][2] + lhs[0][1] * rhs[1][2] + lhs[0][2] * rhs[2][2] + lhs[0][3] * rhs[3][2],
            lhs[0][0] * rhs[0][3] + lhs[0][1] * rhs[1][3] + lhs[0][2] * rhs[2][3] + lhs[0][3] * rhs[3][3],

            lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0] + lhs[1][2] * rhs[2][0] + lhs[1][3] * rhs[3][0],
            lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1] + lhs[1][2] * rhs[2][1] + lhs[1][3] * rhs[3][1],
            lhs[1][0] * rhs[0][2] + lhs[1][1] * rhs[1][2] + lhs[1][2] * rhs[2][2] + lhs[1][3] * rhs[3][2],
            lhs[1][0] * rhs[0][3] + lhs[1][1] * rhs[1][3] + lhs[1][2] * rhs[2][3] + lhs[1][3] * rhs[3][3],

            lhs[2][0] * rhs[0][0] + lhs[2][1] * rhs[1][0] + lhs[2][2] * rhs[2][0] + lhs[2][3] * rhs[3][0],
            lhs[2][0] * rhs[0][1] + lhs[2][1] * rhs[1][1] + lhs[2][2] * rhs[2][1] + lhs[2][3] * rhs[3][1],
            lhs[2][0] * rhs[0][2] + lhs[2][1] * rhs[1][2] + lhs[2][2] * rhs[2][2] + lhs[2][3] * rhs[3][2],
            lhs[2][0] * rhs[0][3] + lhs[2][1] * rhs[1][3] + lhs[2][2] * rhs[2][3] + lhs[2][3] * rhs[3][3],

            lhs[3][0] * rhs[0][0] + lhs[3][1] * rhs[1][0] + lhs[3][2] * rhs[2][0] + lhs[3][3] * rhs[3][0],
            lhs[3][0] * rhs[0][1] + lhs[3][1] * rhs[1][1] + lhs[3][2] * rhs[2][1] + lhs[3][3] * rhs[3][1],
            lhs[3][0] * rhs[0][2] + lhs[3][1] * rhs[1][2] + lhs[3][2] * rhs[2][2] + lhs[3][3] * rhs[3][2],
            lhs[3][0] * rhs[0][3] + lhs[3][1] * rhs[1][3] + lhs[3][2] * rhs[2][3] + lhs[3][3] * rhs[3][3]);
    }

    template <class T>
    void multiply(matrix<T, 4, 4>* out, const matrix<T, 4, 4>& lhs, const matrix<T, 4, 4>& rhs)
    {

        (*out)[0][0] = lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0] + lhs[0][2] * rhs[2][0] + lhs[0][3] * rhs[3][0];
        (*out)[0][1] = lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1] + lhs[0][2] * rhs[2][1] + lhs[0][3] * rhs[3][1];
        (*out)[0][2] = lhs[0][0] * rhs[0][2] + lhs[0][1] * rhs[1][2] + lhs[0][2] * rhs[2][2] + lhs[0][3] * rhs[3][2];
        (*out)[0][3] = lhs[0][0] * rhs[0][3] + lhs[0][1] * rhs[1][3] + lhs[0][2] * rhs[2][3] + lhs[0][3] * rhs[3][3];

        (*out)[1][0] = lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0] + lhs[1][2] * rhs[2][0] + lhs[1][3] * rhs[3][0];
        (*out)[1][1] = lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1] + lhs[1][2] * rhs[2][1] + lhs[1][3] * rhs[3][1];
        (*out)[1][2] = lhs[1][0] * rhs[0][2] + lhs[1][1] * rhs[1][2] + lhs[1][2] * rhs[2][2] + lhs[1][3] * rhs[3][2];
        (*out)[1][3] = lhs[1][0] * rhs[0][3] + lhs[1][1] * rhs[1][3] + lhs[1][2] * rhs[2][3] + lhs[1][3] * rhs[3][3];

        (*out)[2][0] = lhs[2][0] * rhs[0][0] + lhs[2][1] * rhs[1][0] + lhs[2][2] * rhs[2][0] + lhs[2][3] * rhs[3][0];
        (*out)[2][1] = lhs[2][0] * rhs[0][1] + lhs[2][1] * rhs[1][1] + lhs[2][2] * rhs[2][1] + lhs[2][3] * rhs[3][1];
        (*out)[2][2] = lhs[2][0] * rhs[0][2] + lhs[2][1] * rhs[1][2] + lhs[2][2] * rhs[2][2] + lhs[2][3] * rhs[3][2];
        (*out)[2][3] = lhs[2][0] * rhs[0][3] + lhs[2][1] * rhs[1][3] + lhs[2][2] * rhs[2][3] + lhs[2][3] * rhs[3][3];

        (*out)[3][0] = lhs[3][0] * rhs[0][0] + lhs[3][1] * rhs[1][0] + lhs[3][2] * rhs[2][0] + lhs[3][3] * rhs[3][0];
        (*out)[3][1] = lhs[3][0] * rhs[0][1] + lhs[3][1] * rhs[1][1] + lhs[3][2] * rhs[2][1] + lhs[3][3] * rhs[3][1];
        (*out)[3][2] = lhs[3][0] * rhs[0][2] + lhs[3][1] * rhs[1][2] + lhs[3][2] * rhs[2][2] + lhs[3][3] * rhs[3][2];
        (*out)[3][3] = lhs[3][0] * rhs[0][3] + lhs[3][1] * rhs[1][3] + lhs[3][2] * rhs[2][3] + lhs[3][3] * rhs[3][3];

        return;
    }

    //-----------------------------------------------
    // utilities
    template <class T>
    inline matrix<T, 4, 4> transpose(const matrix<T, 4, 4>& rhs)
    {
        return matrix<T, 4, 4>(
            rhs[0][0], rhs[1][0], rhs[2][0], rhs[3][0],
            rhs[0][1], rhs[1][1], rhs[2][1], rhs[3][1],
            rhs[0][2], rhs[1][2], rhs[2][2], rhs[3][2],
            rhs[0][3], rhs[1][3], rhs[2][3], rhs[3][3]);
    }

    template <class T>
    inline T det(const matrix<T, 4, 4>& rhs)
    {
        return rhs.det();
    }

#ifdef minor
//#error ?
#else
    template <class T> //minor
    matrix<T, 4, 4> minor(const matrix<T, 4, 4>& rhs)
    {

        return matrix<T, 4, 4>(
            det33(rhs[1][1], rhs[1][2], rhs[1][3], rhs[2][1], rhs[2][2], rhs[2][3], rhs[3][1], rhs[3][2], rhs[3][3]),
            -det33(rhs[0][1], rhs[0][2], rhs[0][3], rhs[2][1], rhs[2][2], rhs[2][3], rhs[3][1], rhs[3][2], rhs[3][3]),
            det33(rhs[0][1], rhs[0][2], rhs[0][3], rhs[1][1], rhs[1][2], rhs[1][3], rhs[3][1], rhs[3][2], rhs[3][3]),
            -det33(rhs[0][1], rhs[0][2], rhs[0][3], rhs[1][1], rhs[1][2], rhs[1][3], rhs[2][1], rhs[2][2], rhs[2][3]),

            -det33(rhs[1][0], rhs[1][2], rhs[1][3], rhs[2][0], rhs[2][2], rhs[2][3], rhs[3][0], rhs[3][2], rhs[3][3]),
            det33(rhs[0][0], rhs[0][2], rhs[0][3], rhs[2][0], rhs[2][2], rhs[2][3], rhs[3][0], rhs[3][2], rhs[3][3]),
            -det33(rhs[0][0], rhs[0][2], rhs[0][3], rhs[1][0], rhs[1][2], rhs[1][3], rhs[3][0], rhs[3][2], rhs[3][3]),
            det33(rhs[0][0], rhs[0][2], rhs[0][3], rhs[1][0], rhs[1][2], rhs[1][3], rhs[2][0], rhs[2][2], rhs[2][3]),

            det33(rhs[1][0], rhs[1][1], rhs[1][3], rhs[2][0], rhs[2][1], rhs[2][3], rhs[3][0], rhs[3][1], rhs[3][3]),
            -det33(rhs[0][0], rhs[0][1], rhs[0][3], rhs[2][0], rhs[2][1], rhs[2][3], rhs[3][0], rhs[3][1], rhs[3][3]),
            det33(rhs[0][0], rhs[0][1], rhs[0][3], rhs[1][0], rhs[1][1], rhs[1][3], rhs[3][0], rhs[3][1], rhs[3][3]),
            -det33(rhs[0][0], rhs[0][1], rhs[0][3], rhs[1][0], rhs[1][1], rhs[1][3], rhs[2][0], rhs[2][1], rhs[2][3]),

            -det33(rhs[1][0], rhs[1][1], rhs[1][2], rhs[2][0], rhs[2][1], rhs[2][2], rhs[3][0], rhs[3][1], rhs[3][2]),
            det33(rhs[0][0], rhs[0][1], rhs[0][2], rhs[2][0], rhs[2][1], rhs[2][2], rhs[3][0], rhs[3][1], rhs[3][2]),
            -det33(rhs[0][0], rhs[0][1], rhs[0][2], rhs[1][0], rhs[1][1], rhs[1][2], rhs[3][0], rhs[3][1], rhs[3][2]),
            det33(rhs[0][0], rhs[0][1], rhs[0][2], rhs[1][0], rhs[1][1], rhs[1][2], rhs[2][0], rhs[2][1], rhs[2][2]));
    }
#endif

    template <class T> //inverter
    matrix<T, 4, 4> invert(const matrix<T, 4, 4>& rhs)
    {

#if 1
        T det = rhs.det();

        assert(det != T());

        det = T(1) / det;

        return matrix<T, 4, 4>(
            det33(rhs[1][1], rhs[1][2], rhs[1][3], rhs[2][1], rhs[2][2], rhs[2][3], rhs[3][1], rhs[3][2], rhs[3][3]) * det,
            -det33(rhs[0][1], rhs[0][2], rhs[0][3], rhs[2][1], rhs[2][2], rhs[2][3], rhs[3][1], rhs[3][2], rhs[3][3]) * det,
            det33(rhs[0][1], rhs[0][2], rhs[0][3], rhs[1][1], rhs[1][2], rhs[1][3], rhs[3][1], rhs[3][2], rhs[3][3]) * det,
            -det33(rhs[0][1], rhs[0][2], rhs[0][3], rhs[1][1], rhs[1][2], rhs[1][3], rhs[2][1], rhs[2][2], rhs[2][3]) * det,

            -det33(rhs[1][0], rhs[1][2], rhs[1][3], rhs[2][0], rhs[2][2], rhs[2][3], rhs[3][0], rhs[3][2], rhs[3][3]) * det,
            det33(rhs[0][0], rhs[0][2], rhs[0][3], rhs[2][0], rhs[2][2], rhs[2][3], rhs[3][0], rhs[3][2], rhs[3][3]) * det,
            -det33(rhs[0][0], rhs[0][2], rhs[0][3], rhs[1][0], rhs[1][2], rhs[1][3], rhs[3][0], rhs[3][2], rhs[3][3]) * det,
            det33(rhs[0][0], rhs[0][2], rhs[0][3], rhs[1][0], rhs[1][2], rhs[1][3], rhs[2][0], rhs[2][2], rhs[2][3]) * det,

            det33(rhs[1][0], rhs[1][1], rhs[1][3], rhs[2][0], rhs[2][1], rhs[2][3], rhs[3][0], rhs[3][1], rhs[3][3]) * det,
            -det33(rhs[0][0], rhs[0][1], rhs[0][3], rhs[2][0], rhs[2][1], rhs[2][3], rhs[3][0], rhs[3][1], rhs[3][3]) * det,
            det33(rhs[0][0], rhs[0][1], rhs[0][3], rhs[1][0], rhs[1][1], rhs[1][3], rhs[3][0], rhs[3][1], rhs[3][3]) * det,
            -det33(rhs[0][0], rhs[0][1], rhs[0][3], rhs[1][0], rhs[1][1], rhs[1][3], rhs[2][0], rhs[2][1], rhs[2][3]) * det,

            -det33(rhs[1][0], rhs[1][1], rhs[1][2], rhs[2][0], rhs[2][1], rhs[2][2], rhs[3][0], rhs[3][1], rhs[3][2]) * det,
            det33(rhs[0][0], rhs[0][1], rhs[0][2], rhs[2][0], rhs[2][1], rhs[2][2], rhs[3][0], rhs[3][1], rhs[3][2]) * det,
            -det33(rhs[0][0], rhs[0][1], rhs[0][2], rhs[1][0], rhs[1][1], rhs[1][2], rhs[3][0], rhs[3][1], rhs[3][2]) * det,
            det33(rhs[0][0], rhs[0][1], rhs[0][2], rhs[1][0], rhs[1][1], rhs[1][2], rhs[2][0], rhs[2][1], rhs[2][2]) * det);
#else
        matrix<T, 4, 4> t;
        invert(&t, rhs);
        return t;
#endif
    }

    template <class T> //inverter
    inline matrix<T, 4, 4> operator~(const matrix<T, 4, 4>& rhs)
    {
        return invert(rhs);
    }

    template <class T> //inverter
    bool invert_rot_trans(matrix<T, 4, 4>* out, const matrix<T, 4, 4>& rhs)
    {

        //

        //

        (*out)[0][0] = rhs[0][0];
        (*out)[1][0] = rhs[0][1];
        (*out)[2][0] = rhs[0][2];
        (*out)[3][0] = rhs[3][0];

        (*out)[0][1] = rhs[1][0];
        (*out)[1][1] = rhs[1][1];
        (*out)[2][1] = rhs[1][2];
        (*out)[3][1] = rhs[3][1];

        (*out)[0][2] = rhs[2][0];
        (*out)[1][2] = rhs[2][1];
        (*out)[2][2] = rhs[2][2];
        (*out)[3][2] = rhs[3][2];

        (*out)[0][3] = -(rhs[0][0] * rhs[0][3] + rhs[1][0] * rhs[1][3] + rhs[2][0] * rhs[2][3]); //rhs[3][0];
        (*out)[1][3] = -(rhs[0][1] * rhs[0][3] + rhs[1][1] * rhs[1][3] + rhs[2][1] * rhs[2][3]); //rhs[3][1];
        (*out)[2][3] = -(rhs[0][2] * rhs[0][3] + rhs[1][2] * rhs[1][3] + rhs[2][2] * rhs[2][3]); //rhs[3][2];
        (*out)[3][3] = rhs[3][3];

        return true;
    }

    //--------------------------------------------------
    //compare

    template <class T>
    inline bool operator==(const matrix<T, 4, 4>& lhs, const matrix<T, 4, 4>& rhs)
    {
#define DECLARE_OP(i, j) (lhs[i][j] == rhs[i][j])
        return DECLARE_OP(0, 0) && DECLARE_OP(0, 1) && DECLARE_OP(0, 2) && DECLARE_OP(0, 3) &&
               DECLARE_OP(1, 0) && DECLARE_OP(1, 1) && DECLARE_OP(1, 2) && DECLARE_OP(1, 3) &&
               DECLARE_OP(2, 0) && DECLARE_OP(2, 1) && DECLARE_OP(2, 2) && DECLARE_OP(2, 3) &&
               DECLARE_OP(3, 0) && DECLARE_OP(3, 1) && DECLARE_OP(3, 2) && DECLARE_OP(3, 3);
#undef DECLARE_OP
    }

    template <class T>
    inline bool operator!=(const matrix<T, 4, 4>& lhs, const matrix<T, 4, 4>& rhs)
    {
        return !(lhs == rhs);
    }

} //end of namespace

//#include "matrix_functions.hpp"

#endif
