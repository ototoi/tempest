#pragma once
#ifndef TM_VECTOR_FUNCTIONS_EX_HPP
#define TM_VECTOR_FUNCTIONS_EX_HPP

/** 
 * @file vector_functions_ex.hpp
 * @brief Defines operators for anoter type vector.
 *
 * $RCSfile: vector_functions_ex.hpp,v $
 * $Date: 2005/03/13 09:00:42 $
 * $Author: Toru Matsuoka $
 */

#include "operation_result.hpp"

namespace tempest
{

#define DECLARE_OP(OP)                                                                                                       \
    template <class T, class X, std::size_t Sz>                                                                              \
    inline vector<typename operation_result<T, X>::type, Sz> operator OP(const vector<T, Sz>& lhs, const vector<X, Sz>& rhs) \
    {                                                                                                                        \
        typedef vector<typename operation_result<T, X>::type, Sz> vector_type;                                               \
        return vector_type(lhs) OP## = vector_type(rhs);                                                                     \
    }

    DECLARE_OP(+)
    DECLARE_OP(-)
    DECLARE_OP(*)
    DECLARE_OP(/ )

#undef DECLARE_OP

    template <class T, class X, std::size_t Sz>
    inline vector<typename operation_result<T, X>::type, Sz> operator*(const vector<T, Sz>& lhs, X rhs)
    {
        typedef typename operation_result<T, X>::type r_type;
        typedef vector<r_type, Sz> vector_type;
        return vector_type(lhs) *= r_type(rhs);
    }

    template <class T, class X, std::size_t Sz>
    inline vector<typename operation_result<T, X>::type, Sz> operator*(T lhs, const vector<X, Sz>& rhs)
    {
        typedef typename operation_result<T, X>::type r_type;
        typedef vector<r_type, Sz> vector_type;
        return vector_type(rhs) *= r_type(lhs);
    }

    template <class T, class X, std::size_t Sz>
    inline vector<typename operation_result<T, X>::type, Sz> operator/(const vector<T, Sz>& lhs, X rhs)
    {
        typedef typename operation_result<T, X>::type r_type;
        typedef vector<r_type, Sz> vector_type;
        return vector_type(lhs) / r_type(rhs);
    }
}

#endif
