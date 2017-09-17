#pragma once
#ifndef TM_TRANSFORM_HPP
#define TM_TRANSFORM_HPP

/** 
 * @file transform.hpp
 * @brief Defines functions for vector VS. matrix.
 *
 * $RCSfile: transform.hpp,v $
 * $Date: 2005/05/09 12:24:37 $
 * $Author: Toru Matsuoka $
 */

#include <cstddef>

#include "vector_base.hpp"
#include "matrix_base.hpp"

#include "vector.hpp"
#include "matrix.hpp"
#include "quaternion.hpp"

namespace tempest
{

    // Vout = v * m;
    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline vector<T, ColumnSz> operator*(const vector<T, RowSz>& lhs, const matrix<T, RowSz, ColumnSz>& rhs)
    {
        vector<T, ColumnSz> tmp;

        for (std::size_t i = 0; i < ColumnSz; ++i)
        {
            T sum = T();
            for (std::size_t j = 0; j < RowSz; j++)
            {
                sum += lhs[j] * rhs[j][i];
            }
            tmp[i] = sum;
        }
        return tmp;
    }

    // Vout = m * v;
    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline vector<T, RowSz> operator*(const matrix<T, RowSz, ColumnSz>& lhs, const vector<T, ColumnSz>& rhs)
    {
        vector<T, RowSz> tmp;
        for (std::size_t i = 0; i < RowSz; ++i)
        {
            T sum = T();
            for (std::size_t j = 0; j < ColumnSz; j++)
            {
                sum += lhs[i][j] * rhs[j];
            }
            tmp[i] = sum;
        }
        return tmp;
    }
    // v *= m;
    //(v = v * m)
    //template<class T>

    // transform(&out,m,v);
    //(v = m * v)
    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline void transform(
        vector<T, RowSz>* out,
        const matrix<T, RowSz, ColumnSz>& lhs,
        const vector<T, ColumnSz>& rhs)
    {
        vector<T, ColumnSz> tmp;
        for (std::size_t i = 0; i < RowSz; ++i)
        {
            T sum = T();
            for (std::size_t j = 0; j < ColumnSz; j++)
            {
                sum += lhs[i][j] * rhs[j];
            }
            tmp[i] = sum;
        }

        (*out) = tmp;
    }

    //-----------------------------------------------
    //vector<T,2> VS matrix<T,2>

    // Vout = v * m;
    template <class T>
    inline vector<T, 2> operator*(const vector<T, 2>& lhs, const matrix<T, 2>& rhs)
    {
        return vector<T, 2>(lhs[0] * rhs[0][0] + lhs[1] * rhs[1][0], lhs[0] * rhs[0][1] + lhs[1] * rhs[1][1]);
    }
    // Vout = m * v;
    template <class T>
    inline vector<T, 2> operator*(const matrix<T, 2>& lhs, const vector<T, 2>& rhs)
    {
        return vector<T, 2>(lhs[0][0] * rhs[0] + lhs[0][1] * rhs[1], lhs[1][0] * rhs[0] + lhs[1][1] * rhs[1]);
    }
    // v *= m;
    //(v = v * m)
    template <class T>
    inline vector<T, 2>& operator*=(vector<T, 2>& lhs, const matrix<T, 2>& rhs)
    {
        lhs = vector<T, 2>(lhs) * rhs;
        return lhs;
    }
    // transform(m,v);
    //(v = m * v)
    template <class T>
    inline void transform(
        vector<T, 2>* out,
        const matrix<T, 2>& lhs,
        const vector<T, 2>& rhs)
    {
        (*out) = vector<T, 2>(
            lhs[0][0] * rhs[0] + lhs[0][1] * rhs[1],
            lhs[1][0] * rhs[0] + lhs[1][1] * rhs[1]);
    }

    //-----------------------------------------------

    //-----------------------------------------------

    /*
    vector<T,2> VS matrix<T,3>
    */
    // Vout = v * m;
    template <class T>
    inline vector<T, 2> operator*(const vector<T, 2>& lhs, const matrix<T, 3>& rhs)
    {
        return vector<T, 2>(
            lhs[0] * rhs[0][0] + lhs[1] * rhs[1][0] + rhs[2][0],
            lhs[0] * rhs[0][1] + lhs[1] * rhs[1][1] + rhs[2][1]);
    }
    // Vout = m * v;
    template <class T>
    inline vector<T, 2> operator*(const matrix<T, 3>& lhs, const vector<T, 2>& rhs)
    {
        return vector<T, 2>(
            lhs[0][0] * rhs[0] + lhs[0][1] * rhs[1] + lhs[0][2],
            lhs[1][0] * rhs[0] + lhs[1][1] * rhs[1] + lhs[1][2]);
    }
    // v *= m;
    //(v = v * m)
    template <class T>
    inline vector<T, 2>& operator*=(vector<T, 2>& lhs, const matrix<T, 3>& rhs)
    {
        lhs = vector<T, 2>(lhs) * rhs;
        return lhs;
    }
    // transform(m,v);
    //(v = m * v)
    template <class T>
    inline void transform(
        vector<T, 2>* out,
        const matrix<T, 3>& lhs,
        const vector<T, 2>& rhs)
    {
        (*out) = vector<T, 2>(
            lhs[0][0] * rhs[0] + lhs[0][1] * rhs[1] + lhs[0][2],
            lhs[1][0] * rhs[0] + lhs[1][1] * rhs[1] + lhs[1][2]);
    }

    //-----------------------------------------------

    /*
    vector<T,3> VS matrix<T,3>
    */
    // Vout = v * m;
    template <class T>
    inline vector<T, 3> operator*(const vector<T, 3>& lhs, const matrix<T, 3>& rhs)
    {
        return vector<T, 3>(
            lhs[0] * rhs[0][0] + lhs[1] * rhs[1][0] + lhs[2] * rhs[2][0],
            lhs[0] * rhs[0][1] + lhs[1] * rhs[1][1] + lhs[2] * rhs[2][1],
            lhs[0] * rhs[0][2] + lhs[1] * rhs[1][2] + lhs[2] * rhs[2][2]);
    }
    // Vout = m * v;
    template <class T>
    inline vector<T, 3> operator*(const matrix<T, 3>& lhs, const vector<T, 3>& rhs)
    {
        return vector<T, 3>(
            lhs[0][0] * rhs[0] + lhs[0][1] * rhs[1] + lhs[0][2] * rhs[2],
            lhs[1][0] * rhs[0] + lhs[1][1] * rhs[1] + lhs[1][2] * rhs[2],
            lhs[2][0] * rhs[0] + lhs[2][1] * rhs[1] + lhs[2][2] * rhs[2]);
    }
    // v *= m;
    //(v = v * m)
    template <class T>
    inline vector<T, 3>& operator*=(vector<T, 3>& lhs, const matrix<T, 3>& rhs)
    {
        lhs = vector<T, 3>(lhs) * rhs;
        return lhs;
    }
    // transform(m,v);
    //(v = m * v)
    template <class T>
    inline void transform(
        vector<T, 3>* out,
        const matrix<T, 3>& lhs,
        const vector<T, 3>& rhs)
    {
        (*out) = vector<T, 3>(
            lhs[0][0] * rhs[0] + lhs[0][1] * rhs[1] + lhs[0][2] * rhs[2],
            lhs[1][0] * rhs[0] + lhs[1][1] * rhs[1] + lhs[1][2] * rhs[2],
            lhs[2][0] * rhs[0] + lhs[2][1] * rhs[1] + lhs[2][2] * rhs[2]);
    }

    //-----------------------------------------------
    //-----------------------------------------------

    /*
    vector<T,4> VS matrix<T,4>
    */

    // Vout = v * m;
    template <class T>
    inline vector<T, 4> operator*(const vector<T, 4>& lhs, const matrix<T, 4>& rhs)
    {
        return vector<T, 4>(
            lhs[0] * rhs[0][0] + lhs[1] * rhs[1][0] + lhs[2] * rhs[2][0] + lhs[3] * rhs[3][0],
            lhs[0] * rhs[0][1] + lhs[1] * rhs[1][1] + lhs[2] * rhs[2][1] + lhs[3] * rhs[3][1],
            lhs[0] * rhs[0][2] + lhs[1] * rhs[1][2] + lhs[2] * rhs[2][2] + lhs[3] * rhs[3][2],
            lhs[0] * rhs[0][3] + lhs[1] * rhs[1][3] + lhs[2] * rhs[2][3] + lhs[3] * rhs[3][3]);
    }
    // Vout = m * v;
    template <class T>
    inline vector<T, 4> operator*(const matrix<T, 4>& lhs, const vector<T, 4>& rhs)
    {
        return vector<T, 4>(
            lhs[0][0] * rhs[0] + lhs[0][1] * rhs[1] + lhs[0][2] * rhs[2] + lhs[0][3] * rhs[3],
            lhs[1][0] * rhs[0] + lhs[1][1] * rhs[1] + lhs[1][2] * rhs[2] + lhs[1][3] * rhs[3],
            lhs[2][0] * rhs[0] + lhs[2][1] * rhs[1] + lhs[2][2] * rhs[2] + lhs[2][3] * rhs[3],
            lhs[3][0] * rhs[0] + lhs[3][1] * rhs[1] + lhs[3][2] * rhs[2] + lhs[3][3] * rhs[3]);
    }
    // v *= m;
    //(v = v * m)
    template <class T>
    inline vector<T, 4>& operator*=(vector<T, 4>& lhs, const matrix<T, 4>& rhs)
    {
        lhs = vector<T, 4>(lhs) * rhs;
        return lhs;
    }

    // transform(m,v);
    //(v = m * v)
    template <class T>
    inline void transform(
        vector<T, 4>* out,
        const matrix<T, 4>& lhs,
        const vector<T, 4>& rhs)
    {
        (*out) = vector<T, 4>(
            lhs[0][0] * rhs[0] + lhs[0][1] * rhs[1] + lhs[0][2] * rhs[2] + lhs[0][3] * rhs[3],
            lhs[1][0] * rhs[0] + lhs[1][1] * rhs[1] + lhs[1][2] * rhs[2] + lhs[1][3] * rhs[3],
            lhs[2][0] * rhs[0] + lhs[2][1] * rhs[1] + lhs[2][2] * rhs[2] + lhs[2][3] * rhs[3],
            lhs[3][0] * rhs[0] + lhs[3][1] * rhs[1] + lhs[3][2] * rhs[2] + lhs[3][3] * rhs[3]);
    }

    //-----------------------------------------------
    /*
    vector<T,3> VS matrix<T,4>
    */
    // Vout = v * m;
    template <class T>
    inline vector<T, 3> operator*(const vector<T, 3>& lhs, const matrix<T, 4>& rhs)
    {
        vector<T, 4> tmp = vector<T, 4>(lhs[0], lhs[1], lhs[2], T(1)) * rhs;
        T iR = T(1) / tmp[3];
        return vector<T, 3>(tmp[0] * iR, tmp[1] * iR, tmp[2] * iR);
    }
    // Vout = m * v;
    template <class T>
    inline vector<T, 3> operator*(const matrix<T, 4>& lhs, const vector<T, 3>& rhs)
    {
        vector<T, 4> tmp = lhs * vector<T, 4>(rhs[0], rhs[1], rhs[2], T(1));
        T iR = T(1) / tmp[3];
        return vector<T, 3>(tmp[0] * iR, tmp[1] * iR, tmp[2] * iR);
    }
    // v *= m;
    //(v = v * m)
    template <class T>
    inline vector<T, 3>& operator*=(vector<T, 3>& lhs, const matrix<T, 4>& rhs)
    {
        lhs = vector<T, 3>(lhs) * rhs;
        return lhs;
    }

    // transform(m,v);
    //(v = m * v)
    template <class T>
    inline void transform(
        vector<T, 3>* out,
        const matrix<T, 4>& lhs,
        const vector<T, 3>& rhs)
    {
        (*out) = lhs * rhs;
    }

    template <class T, std::size_t Sz>
    inline void transform(
        vector<T, Sz>* out,
        const quaternion<T>& q,
        const vector<T, Sz>& v)
    {

        //out = q * v * q^(-1)
        quaternion<T> tq(q * quaternion<T>::wxyz(T(0), v[0], v[1], v[2]) * ~q); //~q

        (*out)[0] = tq.x();
        (*out)[1] = tq.y();
        (*out)[2] = tq.z();

        return;
    }

} //end of namespace tempest

#endif
