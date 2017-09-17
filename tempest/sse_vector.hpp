#pragma once
#ifndef TM_SSE_VECTOR_HPP
#define TM_SSE_VECTOR_HPP

/**
 * @file sse_vector.hpp
 * @brief Defines tempest::sse_vector.
 * 
 * $RCSfile: sse_vector.hpp,v $
 * $Date: 2005/03/13 09:00:42 $
 * $Author: Toru Matsuoka $
 */

#include <sstream>
#include <xmmintrin.h>
#include <cstddef> //std::size_t
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <cstdlib>

#include "vector_base.hpp"

//__BORLANDC__

#ifdef __GNU_C__
#ifndef __SSE__
#error "F_CK YOU!"
#endif //__SSE__
#endif //__GNU_C__

#define IOF(m, i) reinterpret_cast<value_type*>(&(m))[i]
#define CIOF(m, i) reinterpret_cast<const value_type*>(&(m))[i]

//m128_f32
namespace tempest
{

    template <class Self, class T, std::size_t Sz>
    class sse_vector_base : public vector_base<Self, T, Sz>
    {
    public:
        typedef T value_type;             ///< type of element
        typedef T& reference;             ///< reference
        typedef const T& const_reference; ///< const reference
        typedef Self this_type;           ///< self

        typedef std::size_t size_type;          ///< type of element size
        typedef std::ptrdiff_t difference_type; ///< difference type of pointer

        typedef T* pointer;             ///< pointer type of each element
        typedef const T* const_pointer; ///< const pointer type of each element

        typedef pointer iterator;             ///< iterator
        typedef const_pointer const_iterator; ///< const iterator
    public:
        static const std::size_t c_size = Sz;
        struct no_aligned_tag
        {
        };

    public:
        //-----------------------------------------------
        //capacity
        size_type size() const { return Sz; }
        size_type max_size() const { return Sz; }
        bool empty() const { return false; }
    };

    template <class T, std::size_t Sz>
    class sse_vector;

    template <>
    class sse_vector<float, 4> : public sse_vector_base<sse_vector<float, 4>, float, 4>
    {
    public:
        //-----------------------------------------------
        //functions for iterator
        iterator begin() { return &(IOF(m, 0)); }
        iterator end() { return &(IOF(m, 4)); }
        const_iterator begin() const { return &(CIOF(m, 0)); }
        const_iterator end() const { return &(CIOF(m, 4)); }
        //

    public:
        sse_vector() {}

        sse_vector(const __m128 rhs) : m(rhs) {}

        sse_vector(const sse_vector& rhs) : m(rhs.m) {}

        explicit sse_vector(float a /**/, float b = 0.0f, float c = 0.0f, float d = 0.0f) : m(_mm_setr_ps(a, b, c, d)) {}

        explicit sse_vector(float* array) : m(_mm_load_ps(array)) {}

        sse_vector(float* array, const no_aligned_tag& no_use) : m(_mm_loadu_ps(array)) {}

        ~sse_vector() {}

        //-----------------------------------------------
        //inserters

        this_type& operator=(const sse_vector& rhs)
        {
            m = rhs.m;
            return *this;
        }

        //-----------------------------------------------
        //operators

        this_type& negate()
        {
            //m = _mm_sub_ps(_mm_setzero_ps(),m);
            m = _mm_mul_ps(_mm_set1_ps(-1.0f), m);
            return *this;
        }

        this_type& rcp()
        {
            __m128 r = _mm_rcp_ps(m);
            m = _mm_mul_ps(r, _mm_sub_ps(_mm_set1_ps(2.0f), _mm_mul_ps(m, r)));
            return *this;
        }

        this_type& operator+=(const this_type& rhs)
        {
            m = _mm_add_ps(m, rhs.m);
            return *this;
        }

        this_type& operator-=(const this_type& rhs)
        {
            m = _mm_sub_ps(m, rhs.m);
            return *this;
        }

        this_type& operator*=(const this_type& rhs)
        {
            //__m128 tm = m;
            m = _mm_mul_ps(m, rhs.m);
            return *this;
        }

        this_type& operator/=(const this_type& rhs)
        {
#if 1
            __m128 r = _mm_rcp_ps(rhs.m);
            m = _mm_mul_ps(m, _mm_mul_ps(r, _mm_sub_ps(_mm_set1_ps(2.0f), _mm_mul_ps(rhs.m, r))));
#else
            m = _mm_div_ps(m, rhs.m);
#endif
            return *this;
        }

        this_type& operator*=(float rhs)
        {
            m = _mm_mul_ps(m, _mm_set1_ps(rhs));
            return *this;
        }

        this_type& operator/=(float rhs)
        {
            __m128 x, r;
#if 1
            x = _mm_set1_ps(rhs);

            r = _mm_rcp_ps(x);
            r = _mm_mul_ps(r, _mm_sub_ps(_mm_set1_ps(2.0f), _mm_mul_ps(x, r))); //R(i+1) = R * (2 - x * R)
#else
            x = _mm_set_ss(rhs);

            r = _mm_rcp_ss(x);
            r = _mm_mul_ss(r, _mm_sub_ss(_mm_set_ss(2.0f), _mm_mul_ss(x, r))); //R(i+1) = R * (2 - x * R)
            r = _mm_shuffle_ps(r, r, 0);                                       //ss->ps
#endif
            m = _mm_mul_ps(m, r);
            //m = _mm_div_ps(m,_mm_load1_ps(&rhs));
            return *this;
        }

        //--------------------------------

        value_type& operator[](size_type i)
        {
            return IOF(m, i);
        }

        value_type operator[](size_type i) const
        {
            return CIOF(m, i);
        }

        //-----------------------------------
        //

        __m128 get() const { return m; }
        void set(__m128 rhs) { m = rhs; }

        float* get_ptr() { return reinterpret_cast<float*>(&m); }

        void get(float* array) const { this->geta(array); }       //array must be alighed !
        void geta(float* array) const { _mm_store_ps(array, m); } //array must be alighed !
        void getu(float* array) const { _mm_storeu_ps(array, m); }

        void set(const float* array) { this->seta(array); }         //array must be alighed !
        void seta(const float* array) { m = (_mm_load_ps(array)); } //array must be alighed !
        void setu(const float* array) { m = (_mm_loadu_ps(array)); }

        void assign(const float* a, const float* b)
        {
            assert(b - a <= static_cast<int>(c_size));
            this->set(a);
        }

        //-----------------------------------
        //
        bool equal(const this_type& rhs) const
        {
            return 0xF == _mm_movemask_ps(_mm_cmpeq_ps(m, rhs.m));
        }

        //-----------------------------------------------
        //utilities
        value_type length() const
        {
            return std::sqrt(sqr_length());
        }
        value_type sqr_length() const
        {
            __m128 r = m;
            r = _mm_mul_ps(r, r);
            __m128 x = _mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 0, 3, 2)); //2,1,4,3 3434
            x = _mm_add_ps(r, x);                                     //6,4,6,4 4668
            r = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 3, 0, 1));        //4,6,4,6 6666
            x = _mm_add_ps(r, x);

            return CIOF(x, 0);
        }

        value_type sum() const
        {
            __m128 r = m;
            __m128 x = _mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 0, 3, 2)); // 2, 1, 4, 3
            x = _mm_add_ps(r, x);                                     // 6, 4, 6, 4
            r = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 3, 0, 1));        // 4, 6, 4, 6
            x = _mm_add_ps(r, x);                                     //10,10,10,10 sum of
            return CIOF(x, 0);
        }

        this_type& normalize()
        {
            //T length = sqr_length();         //||V||^2
            //if (length == 0.0f) return *this;
            __m128 r = m;
            r = _mm_mul_ps(r, r);                                     // 4, 3, 2, 1 //x*x
            __m128 x = _mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 0, 3, 2)); // 2, 1, 4, 3
            x = _mm_add_ps(r, x);                                     // 6, 4, 6, 4
            r = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 3, 0, 1));        // 4, 6, 4, 6
            x = _mm_add_ps(r, x);                                     //10,10,10,10 //sum of x*x
                                                                      //if(_mm_comieq_ss(x,_mm_setzero_ps()))return *this;        // if(x == 0)return ;

//------------------------
//
#if 1
            //all ps
            r = _mm_rsqrt_ps(x);
            r = _mm_mul_ps(r, _mm_sub_ps(_mm_set1_ps(1.5f), _mm_mul_ps(_mm_set1_ps(0.5f), _mm_mul_ps(x, _mm_mul_ps(r, r))))); //r = r*(1.5f - 0.5f  * x*r*r)
#else
            //ss
            r = _mm_rsqrt_ss(x);
            r = _mm_mul_ss(r, _mm_sub_ss(_mm_set_ss(1.5f), _mm_mul_ss(_mm_set_ss(0.5f), _mm_mul_ss(x, _mm_mul_ss(r, r))))); //r = r*(1.5f - 0.5f  * x*r*r)
            r = _mm_shuffle_ps(r, r, 0); //ss->ps
#endif
            m = _mm_mul_ps(m, r);

            return *this;
        }

    public:
        static inline __m128 cross_inner(__m128 vec0, __m128 vec1)
        {
            __m128 tmp0, tmp1, tmp2, tmp3;
            tmp0 = _mm_shuffle_ps(vec0, vec0, _MM_SHUFFLE(3, 0, 2, 1));
            tmp1 = _mm_shuffle_ps(vec1, vec1, _MM_SHUFFLE(3, 1, 0, 2));
            tmp2 = _mm_shuffle_ps(vec0, vec0, _MM_SHUFFLE(3, 1, 0, 2));
            tmp3 = _mm_shuffle_ps(vec1, vec1, _MM_SHUFFLE(3, 0, 2, 1));
            return _mm_sub_ps(_mm_mul_ps(tmp2, tmp3), _mm_mul_ps(tmp0, tmp1));
        }

    public:
        const char* debug() const { return "tempest::sse_vector<float,4>"; }

    private:
        __m128 m;
    };

    //-------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------

    template <class T, std::size_t Sz>
    inline sse_vector<T, Sz> operator+(const sse_vector<T, Sz>& rhs)
    {
        return rhs;
    }

    template <class T, std::size_t Sz>
    inline sse_vector<T, Sz> operator-(const sse_vector<T, Sz>& rhs)
    {
        return sse_vector<T, Sz>(rhs).negate();
    }

    template <class T, std::size_t Sz>
    inline sse_vector<T, Sz> rcp(const sse_vector<T, Sz>& rhs)
    {
        return sse_vector<T, Sz>(rhs).rcp();
    }

#define DACLARE_OP(OP)                                                                               \
    template <class T, std::size_t Sz>                                                               \
    inline sse_vector<T, Sz> operator OP(const sse_vector<T, Sz>& lhs, const sse_vector<T, Sz>& rhs) \
    {                                                                                                \
        return sse_vector<T, Sz>(lhs) OP## = rhs;                                                    \
    }

    DACLARE_OP(+)
    DACLARE_OP(-)
    DACLARE_OP(*)
    DACLARE_OP(/ )

#undef DACLARE_OP

    template <class T, std::size_t Sz>
    inline sse_vector<T, Sz> operator*(const sse_vector<T, Sz>& lhs, T rhs)
    {
        return sse_vector<T, Sz>(lhs) *= rhs;
    }

    template <class T, std::size_t Sz>
    inline sse_vector<T, Sz> operator*(T rhs, const sse_vector<T, Sz>& lhs)
    {
        return sse_vector<T, Sz>(lhs) *= rhs;
    }

    template <class T, std::size_t Sz>
    inline sse_vector<T, Sz> operator/(const sse_vector<T, Sz>& lhs, T rhs)
    {
        return sse_vector<T, Sz>(lhs) /= rhs;
    }

    //-----------------------------------------------
    //utility functions
    template <class T, std::size_t Sz>
    inline sse_vector<T, Sz> normalize(const sse_vector<T, Sz>& rhs)
    {
        return sse_vector<T, Sz>(rhs).normalize();
    }

    template <class T, std::size_t Sz>
    inline T length(const sse_vector<T, Sz>& rhs)
    {
        return rhs.length();
    }

    template <class T, std::size_t Sz>
    inline T sqr_length(const sse_vector<T, Sz>& rhs)
    {
        return rhs.sqr_length();
    }

    template <class T, std::size_t Sz>
    inline T sum(const sse_vector<T, Sz>& rhs)
    {
        return rhs.sum();
    }

    template <class T, std::size_t Sz>
    inline T dot(const sse_vector<T, Sz>& lhs, const sse_vector<T, Sz>& rhs)
    {
        return (lhs * rhs).sum();
    }

    template <class T, std::size_t Sz>
    inline sse_vector<T, Sz> cross(const sse_vector<T, Sz>& lhs, const sse_vector<T, Sz>& rhs)
    {
        return sse_vector<T, Sz>(sse_vector<T, Sz>::cross_inner(lhs.get(), rhs.get()));
    }

    //-----------------------------------------------
    //compare functions
    template <class T, std::size_t Sz>
    inline bool operator==(const sse_vector<T, Sz>& lhs, const sse_vector<T, Sz>& rhs)
    {
        return lhs.equal(rhs);
    }

    template <class T, std::size_t Sz>
    inline bool operator!=(const sse_vector<T, Sz>& lhs, const sse_vector<T, Sz>& rhs)
    {
        return !(lhs == rhs);
    }

    /** 
 *    ostream << 
 */
    template <typename T, std::size_t Sz, typename CharT, class Traits>
    std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const sse_vector<T, Sz>& rhs)
    {

        std::basic_ostringstream<CharT, Traits> s;
        s.flags(os.flags());
        s.imbue(os.getloc());
        s.precision(os.precision());
        s << "(";
        for (std::size_t i = 0; i < Sz - 1; ++i)
        {
            s << rhs[i] << ",";
        }
        s << rhs[Sz - 1] << ")";
        return os << s.str();
    }

} //End of namespace

#undef IOF
#undef CIOF

#endif
