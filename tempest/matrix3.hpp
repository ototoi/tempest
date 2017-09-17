#pragma once
#ifndef TM_MATRIX3_HPP
#define TM_MATRIX3_HPP

/**
 * @file matrix3.hpp
 * @brief Defines tempest::matrix<T,3,3>.
 *
 * $RCSfile: matrix3.hpp,v $
 * $Date: 2005/05/22 10:40:00 $
 * $Author: Toru Matsuoka $
 */

//-----------------------------------------------

//#include<cstddef>//std::size_t
//#include<cmath>
//#include<algorithm>//copy
//#include<functional>//plus
#include <cassert>

#include "matrix.hpp"
//#include"matrix_base.hpp"

namespace tempest
{

/**
 *  @class matrix<T,3,3>
 *  @brief matrix class template (Specific 3-dimentional)
 *    
 *  @code
 *  tempest::matrix<double,3> m;    //3-dimentional matrix m
 *  @endcode
 *  @bug    
 *
 */

    template <class T>
    class matrix<T, 3, 3> : public matrix_base<matrix<T, 3, 3>, T, 3, 3>
    {

    public:
        //-----------------------------------------------
        //type defines
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef matrix<T, 3, 3> this_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        typedef T* iterator;
        typedef const T* const_iterator;

    public:
        //-----------------------------------------------
        //size
        static const size_type row_size = 3;
        static const size_type col_size = 3;
        static const size_type c_size = 3 * 3; //9

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
            T _m00, T _m01, T _m02,
            T _m10, T _m11, T _m12,
            T _m20, T _m21, T _m22) : m00(_m00), m01(_m01), m02(_m02),
                                      m10(_m10), m11(_m11), m12(_m12),
                                      m20(_m20), m21(_m21), m22(_m22) {}

        matrix(const matrix<T, 3, 3>& rhs) : m00(rhs.m00), m01(rhs.m01), m02(rhs.m02),
                                             m10(rhs.m10), m11(rhs.m11), m12(rhs.m12),
                                             m20(rhs.m20), m21(rhs.m21), m22(rhs.m22) {}

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
#endif

        template <class X>
        explicit matrix(const matrix<X, 3, 3>& rhs)
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
        explicit matrix(const matrix_base<Self, Type, row_size, col_size>& rhs)
        {
            std::copy(static_cast<const Self&>(rhs).begin(), static_cast<const Self&>(rhs).end(), begin());
        }

        ~matrix() {}

        //swap :-)
        //void swap(this_type& rhs){std::swap(m,rhs.m);}

        //-----------------------------------------------
        //inserters
        template <class X>
        this_type& operator=(const matrix<X, 3, 3>& rhs)
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
            m00 = -m00;
            m01 = -m01;
            m02 = -m02;
            m10 = -m10;
            m11 = -m11;
            m12 = -m12;
            m20 = -m20;
            m21 = -m21;
            m22 = -m22;
            return *this;
        }

#define DECLARE_OP_EQUAL(OP)                     \
    this_type& operator OP(const this_type& rhs) \
    {                                            \
        m00 OP rhs.m00;                          \
        m01 OP rhs.m01;                          \
        m02 OP rhs.m02;                          \
        m10 OP rhs.m10;                          \
        m11 OP rhs.m11;                          \
        m12 OP rhs.m12;                          \
        m20 OP rhs.m20;                          \
        m21 OP rhs.m21;                          \
        m22 OP rhs.m22;                          \
        return *this;                            \
    }

        DECLARE_OP_EQUAL(+= )
        DECLARE_OP_EQUAL(-= )

#undef DECLARE_OP_EQUAL

        //This functon is NOT recommended.
        this_type& operator*=(const this_type& rhs)
        {
            T tp1, tp2;

            tp1 = m00;
            tp2 = m01;

            m00 = tp1 * rhs.m00 + tp2 * rhs.m10 + m02 * rhs.m20;
            m01 = tp1 * rhs.m01 + tp2 * rhs.m11 + m02 * rhs.m21;
            m02 = tp1 * rhs.m02 + tp2 * rhs.m12 + m02 * rhs.m22;

            tp1 = m10;
            tp2 = m11;

            m10 = tp1 * rhs.m00 + tp2 * rhs.m10 + m12 * rhs.m20;
            m11 = tp1 * rhs.m01 + tp2 * rhs.m11 + m12 * rhs.m21;
            m12 = tp1 * rhs.m02 + tp2 * rhs.m12 + m12 * rhs.m22;

            tp1 = m20;
            tp2 = m21;

            m20 = tp1 * rhs.m00 + tp2 * rhs.m10 + m22 * rhs.m20;
            m21 = tp1 * rhs.m01 + tp2 * rhs.m11 + m22 * rhs.m21;
            m22 = tp1 * rhs.m02 + tp2 * rhs.m12 + m22 * rhs.m22;

            return *this;
        }

#define DECLARE_OP_EQUAL_SCALAR(OP)        \
    this_type& operator OP(value_type rhs) \
    {                                      \
        m00 OP rhs;                        \
        m01 OP rhs;                        \
        m02 OP rhs;                        \
        m10 OP rhs;                        \
        m11 OP rhs;                        \
        m12 OP rhs;                        \
        m20 OP rhs;                        \
        m21 OP rhs;                        \
        m22 OP rhs;                        \
        return *this;                      \
    }

        DECLARE_OP_EQUAL_SCALAR(*= )
        DECLARE_OP_EQUAL_SCALAR(/= )

#undef DECLARE_OP_EQUAL_SCALAR

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
            std::swap(m12, m21);

            return *this;
        }

        T det() const
        { //determinant
            //m = 9 p = 5
            return (
                m00 * (m11 * m22 - m12 * m21) +
                m01 * (m12 * m20 - m10 * m22) +
                m02 * (m10 * m21 - m11 * m20));
            /*
        //m = 2*6 = 12 p = 5
        return (
            m00*m11*m22 - m11*m20*m02+
            m01*m12*m20 - m10*m22*m01+
            m02*m10*m21 - m12*m21*m00
        );
        */
        }

        bool is_invertible() const
        {                              // nonsingular
            return (this->det() != 0); //|M| != 0
        }

        const char* debug() const { return "tempest::matrix<T,3,3>"; }

    private:
        union
        {
            struct
            {
                T m00, m01, m02;
                T m10, m11, m12;
                T m20, m21, m22;
            };
            T element[3][3];
        };
    };

    //-----------------------------------------------
    //Not a member!
    //-----------------------------------------------

    //-----------------------------------------------
    //operators

    //
    //+*-/ see matrix.hpp

    template <class T>
    inline T det22(T a, T b, T c, T d)
    {
        return a * d - b * c;
    }

    template <class T>
    matrix<T, 3, 3> operator*(const matrix<T, 3, 3>& lhs, const matrix<T, 3, 3>& rhs)
    {
        return matrix<T, 3, 3>(
            lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0] + lhs[0][2] * rhs[2][0],
            lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1] + lhs[0][2] * rhs[2][1],
            lhs[0][0] * rhs[0][2] + lhs[0][1] * rhs[1][2] + lhs[0][2] * rhs[2][2],

            lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0] + lhs[1][2] * rhs[2][0],
            lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1] + lhs[1][2] * rhs[2][1],
            lhs[1][0] * rhs[0][2] + lhs[1][1] * rhs[1][2] + lhs[1][2] * rhs[2][2],

            lhs[2][0] * rhs[0][0] + lhs[2][1] * rhs[1][0] + lhs[2][2] * rhs[2][0],
            lhs[2][0] * rhs[0][1] + lhs[2][1] * rhs[1][1] + lhs[2][2] * rhs[2][1],
            lhs[2][0] * rhs[0][2] + lhs[2][1] * rhs[1][2] + lhs[2][2] * rhs[2][2]);
    }

    template <class T>
    void multiply(matrix<T, 3, 3>* out, const matrix<T, 3, 3>& lhs, const matrix<T, 3, 3>& rhs)
    {
        (*out)[0][0] = lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0] + lhs[0][2] * rhs[2][0];
        (*out)[0][1] = lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1] + lhs[0][2] * rhs[2][1];
        (*out)[0][2] = lhs[0][0] * rhs[0][2] + lhs[0][1] * rhs[1][2] + lhs[0][2] * rhs[2][2];

        (*out)[1][0] = lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0] + lhs[1][2] * rhs[2][0];
        (*out)[1][1] = lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1] + lhs[1][2] * rhs[2][1];
        (*out)[1][2] = lhs[1][0] * rhs[0][2] + lhs[1][1] * rhs[1][2] + lhs[1][2] * rhs[2][2];

        (*out)[2][0] = lhs[2][0] * rhs[0][0] + lhs[2][1] * rhs[1][0] + lhs[2][2] * rhs[2][0];
        (*out)[2][1] = lhs[2][0] * rhs[0][1] + lhs[2][1] * rhs[1][1] + lhs[2][2] * rhs[2][1];
        (*out)[2][2] = lhs[2][0] * rhs[0][2] + lhs[2][1] * rhs[1][2] + lhs[2][2] * rhs[2][2];
        return;
    }

    //-----------------------------------------------
    // utilities
    template <class T>
    inline matrix<T, 3, 3> transpose(const matrix<T, 3, 3>& rhs)
    {
        return matrix<T, 3, 3>(
            rhs[0][0], rhs[1][0], rhs[2][0],
            rhs[0][1], rhs[1][1], rhs[2][1],
            rhs[0][2], rhs[1][2], rhs[2][2]);
    }

    template <class T>
    inline T det(const matrix<T, 3, 3>& rhs)
    {
        return rhs.det();
    }

    template <class T> //inverter
    inline matrix<T, 3, 3> invert(const matrix<T, 3, 3>& rhs)
    {
        T det = rhs.det();

        assert(det != T());

        det = T(1) / det;

        return matrix<T, 3, 3>(
            det22(rhs[1][1], rhs[1][2], rhs[2][1], rhs[2][2]) * det, -det22(rhs[0][1], rhs[0][2], rhs[2][1], rhs[2][2]) * det, det22(rhs[0][1], rhs[0][2], rhs[1][1], rhs[1][2]) * det,
            -det22(rhs[1][0], rhs[1][2], rhs[2][0], rhs[2][2]) * det, det22(rhs[0][0], rhs[0][2], rhs[2][0], rhs[2][2]) * det, -det22(rhs[0][0], rhs[0][2], rhs[1][0], rhs[1][2]) * det,
            det22(rhs[1][0], rhs[1][1], rhs[2][0], rhs[2][1]) * det, -det22(rhs[0][0], rhs[0][1], rhs[2][0], rhs[2][1]) * det, det22(rhs[0][0], rhs[0][1], rhs[1][0], rhs[1][1]) * det);
    }

    template <class T> //inverter
    inline matrix<T, 3, 3> operator~(const matrix<T, 3, 3>& rhs)
    {
        return invert(rhs);
    }

#ifdef minor
//#error ?
#else
    template <class T> //minor
    inline matrix<T, 3, 3> minor(const matrix<T, 3, 3>& rhs)
    {
        return matrix<T, 3, 3>(
            det22(rhs[1][1], rhs[1][2], rhs[2][1], rhs[2][2]), -det22(rhs[0][1], rhs[0][2], rhs[2][1], rhs[2][2]), det22(rhs[0][1], rhs[0][2], rhs[1][1], rhs[1][2]),
            -det22(rhs[1][0], rhs[1][2], rhs[2][0], rhs[2][2]), det22(rhs[0][0], rhs[0][2], rhs[2][0], rhs[2][2]), -det22(rhs[0][0], rhs[0][2], rhs[1][0], rhs[1][2]),
            det22(rhs[1][0], rhs[1][1], rhs[2][0], rhs[2][1]), -det22(rhs[0][0], rhs[0][1], rhs[2][0], rhs[2][1]), det22(rhs[0][0], rhs[0][1], rhs[1][0], rhs[1][1]));
    }
#endif

    //ex.
    //double det = m.det();
    //if(det != 0){
    //    det = 1.0/det;
    //    im = det*minor(m);
    //}
    //
    //

    //--------------------------------------------------
    //compare

    template <class T>
    inline bool operator==(const matrix<T, 3, 3>& lhs, const matrix<T, 3, 3>& rhs)
    {
        return (lhs[0][0] == rhs[0][0]) && (lhs[0][1] == rhs[0][1]) && (lhs[0][2] == rhs[0][2]) &&
               (lhs[1][0] == rhs[1][0]) && (lhs[1][1] == rhs[1][1]) && (lhs[1][2] == rhs[1][2]) &&
               (lhs[2][0] == rhs[2][0]) && (lhs[2][1] == rhs[2][1]) && (lhs[2][2] == rhs[2][2]);
    }

    template <class T>
    inline bool operator!=(const matrix<T, 3, 3>& lhs, const matrix<T, 3, 3>& rhs)
    {
        return !(lhs == rhs);
    }

} //end of namespace

//#include "matrix_functions.hpp"

#endif
