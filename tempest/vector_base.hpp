#pragma once
#ifndef TM_VECTOR_BASE_HPP
#define TM_VECTOR_BASE_HPP

/** 
 *  @file vector_base.hpp
 *  @brief Defines tempest::vector_base<Self,Type,Sz>.
 *  @auther Toru Matsuoka
 * 
 *  $Id: vector_base.hpp,v 1.4 2005/03/06 14:22:02      Exp $
 */

#include <cstddef>
namespace tempest
{

    /**
    　*  @brief A class template, which is base of vector class.
    　*/
    template <class Self, class Type, std::size_t Sz>
    class vector_base
    {
    };
}

#endif
