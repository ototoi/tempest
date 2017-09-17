#pragma once
#ifndef TM_MATRIX_BASE_HPP
#define TM_MATRIX_BASE_HPP

/**
 *  @file matrix_base.hpp
 *  @brief Defines tempest::matrix_base<Self,Type,RowSz,ColumnSz>.
 *  @author Toru Matsuoka
 *
 *  $Id: matrix_base.hpp,v 1.4 2005/02/26 17:16:18      Exp $
 */

#include <cstddef>

namespace tempest
{

    /**
     *  @brief A class template, which is base of any matrix classes.
     */
    template <class Self, class Type, std::size_t RowSz, std::size_t ColumnSz>
    class matrix_base
    {
    };
}

#endif
