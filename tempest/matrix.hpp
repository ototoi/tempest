#pragma once
#ifndef TM_MATRIX_HPP
#define TM_MATRIX_HPP

/**
 *  @file matrix.hpp
 *  @brief Defines tempest::matrix<T,RowSz,ColumnSz>.
 *  @author Toru Matsuoka
 *  $Id: matrix.hpp,v 1.7 2005/03/06 14:29:21      Exp $
 */

//-----------------------------------------------

//#include<cstddef>//std::size_t
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstring>

#include "matrix_base.hpp"

namespace tempest
{

    namespace detail_m
    {
        template <class T>
        struct is_void
        {
            static const bool value = false;
        };
        template <>
        struct is_void<void>
        {
            static const bool value = true;
        };

        template <class T>
        struct is_pointer
        {
            static const bool value = false;
        };
        template <class T>
        struct is_pointer<T*>
        {
            static const bool value = true;
        };

        template <class T, bool B>
        struct fill_zero_n_imp;

        template <class T>
        struct fill_zero_n_imp<T, true>
        {
            static inline void func(T* first, size_t n)
            {
                using namespace std;
                memset(reinterpret_cast<void*>(first), (int)0, n * sizeof(T));
            }
        };

        template <class T>
        struct fill_zero_n_imp<T, false>
        {
            static inline void func(T* first, size_t n)
            {
                T temp = T();
                while (n--)
                {
                    *first = temp;
                    ++first;
                }
            }
        };
    }

    template <class IT, class Sz>
    inline void fill_zero_n(IT first, Sz n)
    {
        typedef typename std::iterator_traits<IT>::value_type value_type;
        typedef std::numeric_limits<value_type> NL;

        using namespace detail_m;

        static const bool is_optimize =
            (is_pointer<IT>::value) &
            ((NL::is_integer) | (NL::is_iec559) | (is_void<value_type>::value) | (is_pointer<value_type>::value));

        fill_zero_n_imp<value_type, is_optimize>::func(first, n);
    }

    /**
     *  @brief matrix class template.
     *
     *  @code
     *  tempest::matrix<double,5> m;    //5-dimentional matrix m
     *  tempest::matrix<double,10,15> m2;    //10*15 matrix m2
     *  @endcode
     *
     */

    template <class T, std::size_t RowSz, std::size_t ColumnSz = RowSz>
    class matrix : public matrix_base<matrix<T, RowSz, ColumnSz>, T, RowSz, ColumnSz>
    {
    public:
        //-----------------------------------------------
        //type defines
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef matrix<T, RowSz, ColumnSz> this_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        typedef T* pointer;
        typedef const T* const_pointer;

        typedef pointer iterator;
        typedef const_pointer const_iterator;

    public:
        static const size_type row_size = RowSz;
        static const size_type col_size = ColumnSz;
        static const size_type c_size = RowSz * ColumnSz;

    public:
        //-----------------------------------------------
        //functions for iterator
        iterator begin() { return element.begin(); }
        iterator end() { return element.end(); }

        const_iterator begin() const { return element.begin(); }
        const_iterator end() const { return element.end(); }

        //-----------------------------------------------
        //constructors and destructor
        matrix() {}

        matrix(const this_type& rhs) : element(rhs.element) {}

        template <class X>
        explicit matrix(const matrix<X, RowSz, ColumnSz>& rhs)
        {
            std::copy(rhs.begin(), rhs.end(), begin());
        }

        template <class X>
        explicit matrix(const X* rhs)
        {
            std::copy(rhs, rhs + c_size, begin());
        }

        template <class Self, class Type>
        explicit matrix(const matrix_base<Self, Type, RowSz, ColumnSz>& rhs)
        {
            std::copy(static_cast<const Self&>(rhs).begin(), static_cast<const Self&>(rhs).end(), begin());
        }

        ~matrix() {}

        //swap :-)
        //void swap(this_type& rhs){std::swap(m,rhs.m);}

        //-----------------------------------------------
        //inserters
        this_type& operator=(const this_type& rhs)
        {
            element = rhs.element;
            return *this;
        }
        template <class X>
        this_type& operator=(const matrix<X, RowSz, ColumnSz>& rhs)
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
        template <class X>
        this_type& operator+=(const matrix<X, RowSz, ColumnSz>& rhs)
        {
            typename matrix<X, RowSz, ColumnSz>::const_iterator ri = rhs.begin();
            iterator ii = begin();
            iterator ei = end();
            while (ii != ei)
            {
                *ii += static_cast<T>(*ri);
                ++ii;
                ++ri;
            }
            return *this;
        }

        template <class X>
        this_type& operator-=(const matrix<X, RowSz, ColumnSz>& rhs)
        {
            typename matrix<X, RowSz, ColumnSz>::const_iterator ri = rhs.begin();
            iterator ii = begin();
            iterator ei = end();
            while (ii != ei)
            {
                *ii -= static_cast<T>(*ri);
                ++ii;
                ++ri;
            }
            return *this;
        }

        //This functon is NOT recommended.
        template <class X>
        this_type& operator*=(const matrix<X, ColumnSz, ColumnSz>& rhs)
        {
            this_type temp(*this);

            //std::fill_n( (*this).begin(), c_size, T() );
            fill_zero_n((*this).begin(), c_size);

            for (size_type i = 0; i < RowSz; ++i)
            {
                for (size_type k = 0; k < ColumnSz; k++)
                {
                    for (size_type j = 0; j < ColumnSz; ++j)
                    {
                        element[i][j] += static_cast<T>(temp[i][k] * rhs[k][j]);
                    }
                }
            }
            return *this;
        }

        template <class X>
        this_type& operator*=(const X rhs)
        {
            iterator ii = begin();
            iterator ei = end();
            while (ii != ei)
            {
                *ii *= rhs;
                ++ii;
            }
            return *this;
        }
        template <class X>
        this_type& operator/=(const X rhs)
        {
            iterator ii = begin();
            iterator ei = end();
            while (ii != ei)
            {
                *ii /= rhs;
                ++ii;
            }
            return *this;
        }

        T* operator[](size_type i)
        { //pointer what self is const.
            return element.m[i];
        }
        const T* operator[](size_type i) const
        { //pointer what self is const.
            return element.m[i];
        }

        //-----------------------------------------------
        // utilities
        /*
        this_type & transpose(){
            for(std::size_t i = 0;i<Sz;i++){
                for(std::size_t j = i+1;j<Sz;j++){
                    std::swap(element[i][j],element[j][i]);
                }
            }
            return *this;
        }
        */

        /*
        T det() const {//derminant
            return (
                  m00 * det33(m11,m12,m13,m21,m22,m23,m31,m32,m33)
                - m10 * det33(m01,m02,m03,m21,m22,m23,m31,m32,m33)
                + m20 * det33(m01,m02,m03,m11,m12,m13,m31,m32,m33)
                - m30 * det33(m01,m02,m03,m11,m12,m13,m21,m22,m23)
            );
        }


        bool is_invertible() const{// nonsingular
            return (this->det() != 0);//|M| != 0
        }
        */

        const char* debug() const { return "tempest::matrix<T,RowSz,ColumnSz>"; }

    private:
        struct ___m
        {
            T m[RowSz][ColumnSz];

            T* begin() { return &(m[0][0]); }
            T* end() { return &(m[0][0]) + c_size; }
            const T* begin() const { return &(m[0][0]); }
            const T* end() const { return &(m[0][0]) + c_size; }

        } element;
    };

} //end of namespace

#include "matrix_functions.hpp"

#endif
