#pragma once
#ifndef TM_SSE_MATRIX_HPP
#define TM_SSE_MATRIX_HPP

#include <sstream>
#include <xmmintrin.h>
#include <cstddef> //std::size_t
#include <cmath>
#include "matrix_base.hpp"
#include <cstring>

/**
 * @file sse_matrix.hpp
 * @brief Defines tempest::sse_matrix.
 *
 * $RCSfile: sse_matrix.hpp,v $
 * $Date: 2005/03/13 09:00:42 $
 * $Author: Toru Matsuoka $
 */

#ifdef __GNU_C__
#ifndef __SSE__
#error "F_CK YOU!"
#endif //__SSE__
#endif //__GNU_C__

#define IOF(m, i) reinterpret_cast<value_type*>(&(m))[i]
#define CIOF(m, i) reinterpret_cast<const value_type*>(&(m))[i]
#define CCAST(x) (x)

namespace tempest
{

    namespace sse_matrix_util
    {
        inline __m128 SUM_128(__m128 r)
        {                                                             // 4, 3, 2, 1 0123
            __m128 x = _mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 0, 3, 2)); // 2, 1, 4, 3 1032
            x = _mm_add_ps(x, r);                                     // 6, 4, 6, 4 1155
            r = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 3, 0, 1));        // 4, 6, 4, 6 5511
            return _mm_add_ps(r, x);                                  //10,10,10,10 6666
        }
    }

    template <class Self, class T, std::size_t RowSz, std::size_t ColumnSz>
    class sse_matrix_base : public matrix_base<Self, T, RowSz, ColumnSz>
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
        static const std::size_t row_size = RowSz;
        static const std::size_t col_size = ColumnSz;
        static const std::size_t c_size = RowSz * ColumnSz;
        struct no_aligned_tag
        {
        };

    public:
        //-----------------------------------------------
        //capacity
        size_type size() const { return c_size; }
        size_type max_size() const { return c_size; }
        bool empty() const { return false; }
    };

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    class sse_matrix;

    template <>
    class sse_matrix<float, 4, 4> : public sse_matrix_base<sse_matrix<float, 4, 4>, float, 4, 4>
    {
    public:
        //-----------------------------------------------
        //functions for iterator
        iterator begin() { return &(IOF(m[0], 0)); }
        iterator end() { return &(IOF(m[0], 4)); }
        const_iterator begin() const { return &(CIOF(m[3], 0)); }
        const_iterator end() const { return &(CIOF(m[3], 4)); }

    public:
        sse_matrix() {}

        sse_matrix(const sse_matrix& rhs)
        {
#if 1
            m[0] = rhs.m[0];
            m[1] = rhs.m[1];
            m[2] = rhs.m[2];
            m[3] = rhs.m[3];
#else
            std::memcpy(m, rhs.m, sizeof(__m128) * 4);
#endif
        }

        /*
        sse_matrix(__m128 a,__m128 b,__m128 c, __m128 d){
        m[0] = a;
        m[1] = b;
        m[2] = c;
        m[3] = d;
        }
        */

        explicit sse_matrix(const __m128 rhs[4])
        {
            m[0] = rhs[0];
            m[1] = rhs[1];
            m[2] = rhs[2];
            m[3] = rhs[3];
        }

#define INS_MATRIX(NUM) m[NUM] = _mm_setr_ps(_m##NUM##0, _m##NUM##1, _m##NUM##2, _m##NUM##3);
        explicit sse_matrix(
            float _m00, float _m01, float _m02, float _m03,
            float _m10, float _m11, float _m12, float _m13,
            float _m20, float _m21, float _m22, float _m23,
            float _m30, float _m31, float _m32, float _m33)
        {
            INS_MATRIX(0)
            INS_MATRIX(1)
            INS_MATRIX(2)
            INS_MATRIX(3)
        }
#undef INS_MATRIX

        explicit sse_matrix(const float* array)
        {
            m[0] = (_mm_load_ps(CCAST(array + 0)));
            m[1] = (_mm_load_ps(CCAST(array + 4)));
            m[2] = (_mm_load_ps(CCAST(array + 8)));
            m[3] = (_mm_load_ps(CCAST(array + 12)));
        }

        sse_matrix(const float* array, const no_aligned_tag& no_use)
        {
            m[0] = (_mm_loadu_ps(CCAST(array + 0)));
            m[1] = (_mm_loadu_ps(CCAST(array + 4)));
            m[2] = (_mm_loadu_ps(CCAST(array + 8)));
            m[3] = (_mm_loadu_ps(CCAST(array + 12)));
        }

        ~sse_matrix() {}

        //-----------------------------------------------
        //inserters

        this_type& operator=(const sse_matrix& rhs)
        {
            m[0] = rhs.m[0];
            m[1] = rhs.m[1];
            m[2] = rhs.m[2];
            m[3] = rhs.m[3];
            return *this;
        }

        //-----------------------------------------------
        //operators

        this_type& negate()
        {
            /*
            m[0] = _mm_sub_ps(_mm_setzero_ps(),m[0]);
            m[1] = _mm_sub_ps(_mm_setzero_ps(),m[1]);
            m[2] = _mm_sub_ps(_mm_setzero_ps(),m[2]);
            m[3] = _mm_sub_ps(_mm_setzero_ps(),m[3]);
            */
            m[0] = _mm_mul_ps(_mm_set1_ps(-1.0f), m[0]);
            m[1] = _mm_mul_ps(_mm_set1_ps(-1.0f), m[1]);
            m[2] = _mm_mul_ps(_mm_set1_ps(-1.0f), m[2]);
            m[3] = _mm_mul_ps(_mm_set1_ps(-1.0f), m[3]);

            return *this;
        }

        this_type& operator+=(const this_type& rhs)
        {
            m[0] = _mm_add_ps(m[0], rhs.m[0]);
            m[1] = _mm_add_ps(m[1], rhs.m[1]);
            m[2] = _mm_add_ps(m[2], rhs.m[2]);
            m[3] = _mm_add_ps(m[3], rhs.m[3]);
            return *this;
        }

        this_type& operator-=(const this_type& rhs)
        {
            m[0] = _mm_sub_ps(m[0], rhs.m[0]);
            m[1] = _mm_sub_ps(m[1], rhs.m[1]);
            m[2] = _mm_sub_ps(m[2], rhs.m[2]);
            m[3] = _mm_sub_ps(m[3], rhs.m[3]);
            return *this;
        }

        /*
          [a0,a1,a2,a3] [x0,x1,x2,x3]
          [b0,b1,b2,b3]*[y0,y1,y2,y3]
          [c0,c1,c2,c3] [z0,z1,z2,z3]
          [d0,d1,d2,d3] [w0,w1,w2,w3]
          =
          a0*x0 + a1*y0 + a2*z0 + a3*w0   a0*x1 + a1*y1 + a2*z1 + a3*w1
          a0(x0,x1,x2,x3) + a1(y0,y1,y2,y3)
        */

        this_type& operator*=(const this_type& rhs)
        {

//for(int i=0;i<4;i++){
#define DACLARE_MUL_IMP(i)                                                                  \
    m[i] =                                                                                  \
        _mm_add_ps(                                                                         \
            _mm_add_ps(                                                                     \
                _mm_mul_ps(_mm_shuffle_ps(m[i], m[i], _MM_SHUFFLE(0, 0, 0, 0)), rhs.m[0]),  \
                _mm_mul_ps(_mm_shuffle_ps(m[i], m[i], _MM_SHUFFLE(1, 1, 1, 1)), rhs.m[1])), \
            _mm_add_ps(                                                                     \
                _mm_mul_ps(_mm_shuffle_ps(m[i], m[i], _MM_SHUFFLE(2, 2, 2, 2)), rhs.m[2]),  \
                _mm_mul_ps(_mm_shuffle_ps(m[i], m[i], _MM_SHUFFLE(3, 3, 3, 3)), rhs.m[3])));

            DACLARE_MUL_IMP(0)
            DACLARE_MUL_IMP(1)
            DACLARE_MUL_IMP(2)
            DACLARE_MUL_IMP(3)

#undef DACLARE_MUL_IMP
            //}

            return *this;
        }

        this_type& operator*=(float rhs)
        {
            __m128 rr = _mm_set1_ps(rhs);
            m[0] = _mm_mul_ps(m[0], rr);
            m[1] = _mm_mul_ps(m[1], rr);
            m[2] = _mm_mul_ps(m[2], rr);
            m[3] = _mm_mul_ps(m[3], rr);

            return *this;
        }

        this_type& operator/=(float rhs)
        {
            __m128 x, r;
#if 1
            x = _mm_set1_ps(rhs);

            r = _mm_rcp_ps(x);
            r = _mm_mul_ps(r, _mm_sub_ps(_mm_set1_ps(2), _mm_mul_ps(x, r))); //R(i+1) = R * (2 - x * R)
#else
            x = _mm_set_ss(rhs);

            r = _mm_rcp_ss(x);
            r = _mm_mul_ss(r, _mm_sub_ss(_mm_set_ss(2), _mm_mul_ss(x, r))); //R(i+1) = R * (2 - x * R)
            r = _mm_shuffle_ps(r, r, 0); //ss->ps
#endif

            m[0] = _mm_mul_ps(m[0], r);
            m[1] = _mm_mul_ps(m[1], r);
            m[2] = _mm_mul_ps(m[2], r);
            m[3] = _mm_mul_ps(m[3], r);
            //m = _mm_div_ps(m,_mm_load1_ps(&rhs));
            return *this;
        }

        //--------------------------------

        float* operator[](size_type i)
        {
            return reinterpret_cast<float*>(&(m[i]));
        }

        const float* operator[](size_type i) const
        {
            return reinterpret_cast<const float*>(&(m[i]));
        }

        //-----------------------------------
        //

        const __m128& get(int i) const { return m[i]; }
        void set(int i, __m128 rhs) { m[i] = rhs; }

        void get(float* array) const
        { //array must be alighed !
            this->geta(array);
        }
        void geta(float* array) const
        { //array must be alighed !
            _mm_store_ps(array + 0, m[0]);
            _mm_store_ps(array + 4, m[1]);
            _mm_store_ps(array + 8, m[2]);
            _mm_store_ps(array + 12, m[3]);
        }
        void getu(float* array) const
        { //array must be alighed !
            _mm_storeu_ps(array + 0, m[0]);
            _mm_storeu_ps(array + 4, m[1]);
            _mm_storeu_ps(array + 8, m[2]);
            _mm_storeu_ps(array + 12, m[3]);
        }

        void set(const float* array)
        { //array must be alighed !
            this->seta(array);
        }
        void seta(const float* array)
        { //array must be alighed !
            m[0] = (_mm_load_ps(CCAST(array + 0)));
            m[1] = (_mm_load_ps(CCAST(array + 4)));
            m[2] = (_mm_load_ps(CCAST(array + 8)));
            m[3] = (_mm_load_ps(CCAST(array + 12)));
        }
        void setu(const float* array)
        {
            m[0] = (_mm_loadu_ps(CCAST(array + 0)));
            m[1] = (_mm_loadu_ps(CCAST(array + 4)));
            m[2] = (_mm_loadu_ps(CCAST(array + 8)));
            m[3] = (_mm_loadu_ps(CCAST(array + 12)));
        }

        void assign(const float* a, const float* b)
        {
            assert(b - a <= static_cast<int>(c_size));
            this->set(a);
        }
        //-----------------------------------
        //
        static bool equal(__m128 a, __m128 b)
        {
            return 0xF == _mm_movemask_ps(_mm_cmpeq_ps(a, b));
        }

        bool equal(const this_type& rhs) const
        {
            return equal(m[0], rhs.m[0]) && equal(m[1], rhs.m[1]) && equal(m[2], rhs.m[2]) && equal(m[3], rhs.m[3]);
        }

        //-----------------------------------------------
        // utilities
        this_type& transpose()
        {
            _MM_TRANSPOSE4_PS(m[0], m[1], m[2], m[3]);
            return *this;
        }

        float det() const
        { //derminant
            using namespace sse_matrix_util;

            //a/b/c
            __m128 a[3], b[3];
            __m128 x;

            //m ...3210
            //---------------------------------------
            //0

            a[0] = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(1, 2, 3, 0)); //1
            a[1] = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(2, 3, 1, 0)); //2
            a[2] = _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(3, 1, 2, 0)); //3
            //
            b[0] = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(3, 1, 2, 0)); //1
            b[1] = a[1];                                                //2
            b[2] = _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(1, 2, 3, 0)); //3

            x =
                _mm_mul_ps(
                    _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(0, 0, 0, 0)),
                    _mm_sub_ps(_mm_mul_ps(a[2], _mm_mul_ps(a[1], a[0])), _mm_mul_ps(b[2], _mm_mul_ps(b[1], b[0]))) //(a-b)
                    );
            //+

            //1
            a[0] = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(1, 2, 3, 0)); //0
            //a[1] = _mm_shuffle_ps(m[2],m[2],_MM_SHUFFLE(2,3,1,0));//2
            //a[2] = _mm_shuffle_ps(m[3],m[3],_MM_SHUFFLE(3,1,2,0));//3
            //
            b[0] = _mm_shuffle_ps(m[0], m[0], _MM_SHUFFLE(3, 1, 2, 0)); //0
            //b[1] = _mm_shuffle_ps(m[2],m[2],_MM_SHUFFLE(2,3,1,0));//2
            //b[2] = _mm_shuffle_ps(m[3],m[3],_MM_SHUFFLE(1,2,3,0));//3

            x =
                _mm_sub_ps(
                    x,
                    _mm_mul_ps(
                        _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 0)),
                        _mm_sub_ps(_mm_mul_ps(a[2], _mm_mul_ps(a[1], a[0])), _mm_mul_ps(b[2], _mm_mul_ps(b[1], b[0]))) //(a-b)
                        ));                                                                                            //-

            //2
            //a[0] = _mm_shuffle_ps(m[0],m[0],_MM_SHUFFLE(1,2,3,0));//0
            a[1] = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(2, 3, 1, 0)); //1
            //a[2] = _mm_shuffle_ps(m[3],m[3],_MM_SHUFFLE(3,1,2,0));//3
            //
            //b[0] = _mm_shuffle_ps(m[0],m[0],_MM_SHUFFLE(3,1,2,0));//0
            b[1] = a[1]; //1
            //b[2] = _mm_shuffle_ps(m[3],m[3],_MM_SHUFFLE(1,2,3,0));//3

            x =
                _mm_add_ps(
                    x,
                    _mm_mul_ps(
                        _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 0, 0)),
                        _mm_sub_ps(_mm_mul_ps(a[2], _mm_mul_ps(a[1], a[0])), _mm_mul_ps(b[2], _mm_mul_ps(b[1], b[0]))) //(a-b)
                        ));                                                                                            //+

            //3
            //a[0] = _mm_shuffle_ps(m[0],m[0],_MM_SHUFFLE(1,2,3,0));//0
            //a[1] = _mm_shuffle_ps(m[1],m[1],_MM_SHUFFLE(2,3,1,0));//1
            a[2] = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(3, 1, 2, 0)); //2
            //
            //b[0] = _mm_shuffle_ps(m[0],m[0],_MM_SHUFFLE(3,1,2,0));//0
            //b[1] = _mm_shuffle_ps(m[1],m[1],_MM_SHUFFLE(2,3,1,0));//1
            b[2] = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(1, 2, 3, 0)); //2

            x =
                _mm_sub_ps(
                    x,
                    _mm_mul_ps(
                        _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(0, 0, 0, 0)),
                        _mm_sub_ps(_mm_mul_ps(a[2], _mm_mul_ps(a[1], a[0])), _mm_mul_ps(b[2], _mm_mul_ps(b[1], b[0]))) //(a-b)
                        ));                                                                                            //-

            //---------------------------------------
            IOF(x, 0) = 0.0f; //last is 0 ../ 0123
            x = SUM_128(x);

            return IOF(x, 0);
        }
//...

#define DEC(K_, J_)                                                                                      \
    \
u = _mm_sub_ps(_mm_setzero_ps(), _mm_mul_ps(z[J_], r)); /*  -IOF(z[J_],K_)*IOF(r,K_)  */                 \
    \
z[J_] = _mm_sub_ps(z[J_], _mm_mul_ps(z[K_], _mm_shuffle_ps(z[J_], z[J_], _MM_SHUFFLE(K_, K_, K_, K_)))); \
    \
IOF(z[J_], K_) = IOF(u, K_);

#define DIV_RCP(RET, XXX)                                          \
    \
r = _mm_rcp_ps(XXX);                                               \
    \
r = _mm_mul_ps(r, _mm_sub_ps(_mm_set1_ps(2), _mm_mul_ps(XXX, r))); \
    \
RET = _mm_mul_ps(RET, r);
        this_type inv() const
        {
            __m128 t, u, r, z[4];
            z[0] = get(0);
            z[1] = get(1);
            z[2] = get(2);
            z[3] = get(3);

            //k = 0;
            //tx = IOF(z[0],0);
            t = _mm_shuffle_ps(z[0], z[0], _MM_SHUFFLE(0, 0, 0, 0));
            DIV_RCP(z[0], t)
            //
            IOF(z[0], 0) = IOF(r, 0);
            //------
            DEC(0, 1)
            DEC(0, 2)
            DEC(0, 3)

            //k = 1;
            //tx = IOF(z[1],1);
            t = _mm_shuffle_ps(z[1], z[1], _MM_SHUFFLE(1, 1, 1, 1));
            DIV_RCP(z[1], t)
            IOF(z[1], 1) = IOF(r, 1);
            //------
            DEC(1, 0)
            DEC(1, 2)
            DEC(1, 3)

            //k = 2;
            //tx = IOF(z[2],2);
            t = _mm_shuffle_ps(z[2], z[2], _MM_SHUFFLE(2, 2, 2, 2));
            DIV_RCP(z[2], t)
            IOF(z[2], 2) = IOF(r, 2);
            //------
            DEC(2, 0)
            DEC(2, 1)
            DEC(2, 3)

            //k = 3;
            //tx = IOF(z[3],3);
            t = _mm_shuffle_ps(z[3], z[3], _MM_SHUFFLE(3, 3, 3, 3));
            DIV_RCP(z[3], t)
            IOF(z[3], 3) = IOF(r, 3);
            //------
            DEC(3, 0)
            DEC(3, 1)
            DEC(3, 2)

            return sse_matrix(z);
        }

#undef DEC
#undef DIV_RCP

        bool is_invertible() const
        {                              // nonsingular
            return (this->det() != 0); //|M| != 0
        }

    public:
        const char* debug() const { return "tempest::sse_matrix<float,4,4>"; }

    private: //not privat
        __m128 m[4];
    };

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline sse_matrix<T, RowSz, ColumnSz> operator+(const sse_matrix<T, RowSz, ColumnSz>& rhs)
    {
        return rhs;
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline sse_matrix<T, RowSz, ColumnSz> operator-(const sse_matrix<T, RowSz, ColumnSz>& rhs)
    {
        return sse_matrix<T, RowSz, ColumnSz>(rhs).negate();
    }

#define DACLARE_OP(OP)                                                                                \
    \
template<class T, std::size_t RowSz, std::size_t ColumnSz> \
inline sse_matrix<T, RowSz, ColumnSz>                                                                 \
    operator OP(const sse_matrix<T, RowSz, ColumnSz>& lhs, const sse_matrix<T, RowSz, ColumnSz>& rhs) \
    {                                                                                                 \
        return sse_matrix<T, RowSz, ColumnSz>(lhs) OP## = rhs;                                        \
    \
}

    DACLARE_OP(+)
    DACLARE_OP(-)
    DACLARE_OP(*)

#undef DACLARE_OP

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline sse_matrix<T, RowSz, ColumnSz> operator*(const sse_matrix<T, RowSz, ColumnSz>& lhs, T rhs)
    {
        return sse_matrix<T, RowSz, ColumnSz>(lhs) *= rhs;
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline sse_matrix<T, RowSz, ColumnSz> operator*(T rhs, const sse_matrix<T, RowSz, ColumnSz>& lhs)
    {
        return sse_matrix<T, RowSz, ColumnSz>(lhs) *= rhs;
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline sse_matrix<T, RowSz, ColumnSz> operator/(const sse_matrix<T, RowSz, ColumnSz>& lhs, T rhs)
    {
        return sse_matrix<T, RowSz, ColumnSz>(lhs) /= rhs;
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline sse_matrix<T, RowSz, ColumnSz> operator~(const sse_matrix<T, RowSz, ColumnSz>& rhs)
    {
        return rhs.inv();
    }

    //-----------------------------------------------
    //compare
    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline bool operator==(const sse_matrix<T, RowSz, ColumnSz>& lhs, const sse_matrix<T, RowSz, ColumnSz>& rhs)
    {
        return lhs.equal(rhs);
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline bool operator!=(const sse_matrix<T, RowSz, ColumnSz>& lhs, const sse_matrix<T, RowSz, ColumnSz>& rhs)
    {
        return !(lhs == rhs);
    }

    //-----------------------------------------------
    //output

    /**
     * ostream <<
     */
    template <typename T, std::size_t RowSz, std::size_t ColumnSz, typename _CharT, class _Traits>
    std::basic_ostream<_CharT, _Traits>& operator<<(std::basic_ostream<_CharT, _Traits>& os, const sse_matrix<T, RowSz, ColumnSz>& rhs)
    {

        std::basic_ostringstream<_CharT, _Traits> s;
        s.flags(os.flags());
        s.imbue(os.getloc());
        s.precision(os.precision());
        s << "(";
        for (std::size_t i = 0; i < RowSz; i++)
        {
            s << "(";
            for (std::size_t j = 0; j < ColumnSz; ++j)
            {
                s << rhs[i][j];
                if (j != ColumnSz - 1)
                {
                    s << ",";
                }
            }
            s << ")";
            if (i != RowSz - 1)
            {
                s << ",";
            }
        }
        s << ")";
        return os << s.str();
    }

} //End of namespace

#undef IOF
#undef CIOF

#undef CCAST

#endif
