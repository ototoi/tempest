#pragma once
#ifndef TM_MATRIX_FUNCTIONS_EX_HPP
#define TM_MATRIX_FUNCTIONS_EX_HPP

/** 
 *  @file another_matrix_functions.hpp
 *  @brief Defines operators for anoter type matrix.
 *  @author Toru Matsuoka
 *  $Id: another_matrix_functions.hpp,v 1.5 2005/03/13 09:00:42      Exp $
 */

#include "operation_result.hpp"

namespace tempest
{

#define DECLARE_OP(OP)                                                                   \
    template <class T, class X, std::size_t RowSz, std::size_t ColSz>                    \
    inline matrix<typename operation_result<T, X>::type, RowSz, ColSz>                   \
    operator OP(const matrix<T, RowSz, ColSz>& lhs, const matrix<X, RowSz, ColSz>& rhs)  \
    {                                                                                    \
        typedef matrix<typename operation_result<T, X>::type, RowSz, ColSz> matrix_type; \
        return matrix_type(lhs) OP## = matrix_type(rhs);                                 \
    }

    DECLARE_OP(+)
    DECLARE_OP(-)
    DECLARE_OP(*)
//DECLARE_OP(/)

#undef DECLARE_OP

    template <class T, class X, std::size_t RowSz, std::size_t ColSz>
    inline matrix<typename operation_result<T, X>::type, RowSz, ColSz>
    operator*(const matrix<T, RowSz, ColSz>& lhs, X rhs)
    {
        typedef typename operation_result<T, X>::type r_type;
        typedef matrix<r_type, RowSz, ColSz> matrix_type;
        return matrix_type(lhs) *= r_type(rhs);
    }

    template <class T, class X, std::size_t RowSz, std::size_t ColSz>
    inline matrix<typename operation_result<T, X>::type, RowSz, ColSz>
    operator*(T lhs, const matrix<X, RowSz, ColSz>& rhs)
    {
        typedef typename operation_result<T, X>::type r_type;
        typedef matrix<r_type, RowSz, ColSz> matrix_type;
        return matrix_type(rhs) *= r_type(lhs);
    }

    template <class T, class X, std::size_t RowSz, std::size_t ColSz>
    inline matrix<typename operation_result<T, X>::type, RowSz, ColSz>
    operator/(const matrix<T, RowSz, ColSz>& lhs, X rhs)
    {
        typedef typename operation_result<T, X>::type r_type;
        typedef matrix<r_type, RowSz, ColSz> matrix_type;
        return matrix_type(lhs) /= r_type(rhs);
    }
}

#endif
