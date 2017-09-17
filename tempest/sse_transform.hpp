#pragma once
#ifndef TM_SSE_TRANSFORM_HPP
#define TM_SSE_TRANSFORM_HPP

/** 
 * @file sse_transform.hpp
 * @brief Defines tempest::sse_vector operators.
 *
 * $RCSfile: sse_transform.hpp,v $
 * $Date: 2005/03/06 14:22:02 $
 * $Author: Toru Matsuoka $
 */

#include <cstdlib>
#include "sse_vector.hpp"
#include "sse_matrix.hpp"

namespace tempest
{

    //1x4 * 4x4 -> 1x4
    /**
     *  vout = v*m
     *
     *                          [a,b,c,d]
     *  [X,Y,Z,W] = [x,y,z,w] * [e,f,g,h]  
     *                          [i,j,k,l]
     *                          [m,n,o,p]
     */
    inline sse_vector<float, 4> operator*(const sse_vector<float, 4>& lhs, const sse_matrix<float, 4, 4>& rhs)
    {
        __m128 x = lhs.get();

        //1
        return sse_vector<float, 4>(
            _mm_add_ps(
                _mm_add_ps(
                    _mm_mul_ps(_mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 0, 0)), rhs.get(0)),
                    _mm_mul_ps(_mm_shuffle_ps(x, x, _MM_SHUFFLE(1, 1, 1, 1)), rhs.get(1))),
                _mm_add_ps(
                    _mm_mul_ps(_mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 2, 2, 2)), rhs.get(2)),
                    _mm_mul_ps(_mm_shuffle_ps(x, x, _MM_SHUFFLE(3, 3, 3, 3)), rhs.get(3)))));
    }

    /**
    *  @note 4x4 * 4x1 -> 4x1
    *
    *  vout = m*v
    *
    *  [a,b,c,d]   [x]
    *  [e,f,g,h] * [y]
    *  [i,j,k,l]   [z]
    *  [m,n,o,p]   [w]
    */
    inline sse_vector<float, 4> operator*(const sse_matrix<float, 4, 4>& lhs, const sse_vector<float, 4>& rhs)
    {

        __m128 a = _mm_mul_ps(lhs.get(0), rhs.get());
        __m128 b = _mm_mul_ps(lhs.get(1), rhs.get());
        __m128 c = _mm_mul_ps(lhs.get(2), rhs.get());
        __m128 d = _mm_mul_ps(lhs.get(3), rhs.get());

        _MM_TRANSPOSE4_PS(a, b, c, d);

        return sse_vector<float, 4>(_mm_add_ps(_mm_add_ps(a, b), _mm_add_ps(c, d)));
    }

    inline void transform(
        sse_vector<float, 4>* out,
        const sse_matrix<float, 4, 4>& lhs, const sse_vector<float, 4>& rhs)
    {
        *out = lhs * rhs;
    }
}

#endif
