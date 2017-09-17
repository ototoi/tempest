#pragma once
#ifndef TM_QUATERNION_HPP
#define TM_QUATERNION_HPP

/** 
 *  @file quaternion.hpp
 *  @brief Defines tempest::quaternion<T>.
 *  @author Toru Matsuoka
 * 
 *  $Id: quaternion.hpp,v 1.7 2005/03/17 13:00:40      Exp $
 *
 */

/*
 caution:
 quaternion's inside element 4 vector is Not decideed order "w,x,y,z" or "x,y,z,w".

 example..
 typedef tempest::vector<float,3> VECTOR;
 typedef tempest::quaternion<float> QT;
 QT q1 = QT::wxyz(1,0,0,0);//wxyz 
 QT q2 = QT::xyzw(0,0,0,1);//xyzw
 QT q3 = pivot_angle(VECTOR(1,1,0),30.f*3.1415.f/180.f);
 
 */

#include <cstddef>
#include <cmath>
#include <stdexcept>

#include <complex>

#include <ostream>
#include <sstream>

#include "vector_base.hpp"

namespace tempest
{

    template <class T>
    class quaternion;

    /**
     *  @brief Quaternion class template.
     *  
     *  @code
     *  typedef tempest::vector<float,3> VECTOR;
     *  typedef tempest::quaternion<float> QT;
     *  QT q1 = QT::wxyz(1,0,0,0);//wxyz 
     *  QT q2 = QT::xyzw(0,0,0,1);//xyzw
     *  QT q3 = pivot_angle(VECTOR(1,1,0),30.f*3.1415.f/180.f);// 
     *  @endcode    
     *
     */
    template <class T>
    class quaternion : public vector_base<quaternion<T>, T, 4>
    {
    public:
        //-----------------------------------------------
        //type defines
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef quaternion<T> this_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        typedef T* pointer;
        typedef const T* const_pointer;

        typedef pointer iterator;
        typedef const_pointer const_iterator;

    public:
        static const size_type c_size = 4; //container size

        enum
        {
            W_INDEX, //0
            X_INDEX, //1
            Y_INDEX, //2
            Z_INDEX  //3
        };

        static inline this_type wxyz(value_type w, value_type x, value_type y, value_type z) { return this_type(w, x, y, z); }
        static inline this_type xyzw(value_type x, value_type y, value_type z, value_type w) { return this_type(w, x, y, z); }

    public:
        //-----------------------------------------------
        //functions for iterator
        iterator begin() { return element; }
        iterator end() { return element + c_size; }
        const_iterator begin() const { return element; }
        const_iterator end() const { return element + c_size; }

        //-----------------------------------------------
        //constructors and destructor
        quaternion() {}

        quaternion(value_type a, value_type b, value_type c, value_type d) : m0(a), m1(b), m2(c), m3(d) {}

        quaternion(const quaternion<T>& rhs) : m0(rhs.m0), m1(rhs.m1), m2(rhs.m2), m3(rhs.m3) {}

        template <class X>
        explicit quaternion(const X* rhs) : m0(rhs[0]), m1(rhs[1]), m2(rhs[2]), m3(rhs[3])
        {
        } //X const (&rhs)[c_size]

        template <class X>
        explicit quaternion(const quaternion<X>& rhs) : m0(rhs[0]), m1(rhs[1]), m2(rhs[2]), m3(rhs[3])
        {
        }

        explicit quaternion(
            std::complex<T> const& z0,
            std::complex<T> const& z1 = std::complex<T>()) : m0(z0.real()),
                                                             m1(z0.imag()),
                                                             m2(z1.real()),
                                                             m3(z1.imag()) {} // nothing to do!

        template <class Self, class X>
        explicit quaternion(const vector_base<Self, X, 4>& rhs) : m0(static_cast<T>(static_cast<const Self&>(rhs)[0])),
                                                                  m1(static_cast<T>(static_cast<const Self&>(rhs)[1])),
                                                                  m2(static_cast<T>(static_cast<const Self&>(rhs)[2])),
                                                                  m3(static_cast<T>(static_cast<const Self&>(rhs)[3]))
        {
        }

        //-----------------------------------------------
        //inserters
        this_type& operator=(const this_type& rhs)
        {
            m0 = rhs.m0;
            m1 = rhs.m1;
            m2 = rhs.m2;
            m3 = rhs.m3;
            return *this;
        }

        template <class X>
        this_type& operator=(const quaternion<X>& rhs)
        {
            m0 = rhs[0];
            m1 = rhs[1];
            m2 = rhs[2];
            m3 = rhs[3];
            return *this;
        }

        template <class Self, class X>
        this_type& operator=(const vector_base<Self, X, 4>& rhs)
        {
            m0 = static_cast<T>(static_cast<const Self&>(rhs)[0]);
            m1 = static_cast<T>(static_cast<const Self&>(rhs)[1]);
            m2 = static_cast<T>(static_cast<const Self&>(rhs)[2]);
            m3 = static_cast<T>(static_cast<const Self&>(rhs)[3]);
            return *this;
        }

        template <class IT>
        void assign(IT start, IT end)
        {
            assert(std::distance(start, end) <= c_size); //debug

            std::copy(start, end, begin());
        }

        void assign(size_type num, value_type val)
        {
            //std::fill_n(begin(),(num<c_size)?num:c_size,val);
            size_t sz = (num < c_size) ? num : c_size;
            for (size_t i = 0; i < sz; i++)
                element[i] = val;
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
            m0 = -m0;
            m1 = -m1;
            m2 = -m2;
            m3 = -m3;

            return *this;
        }

// for T
#define DECLARE_OP_EQUAL(OP)                           \
    template <class X>                                 \
    this_type& operator OP(const quaternion<X>& rhs)   \
    {                                                  \
        m0 OP rhs.m0;                                  \
        m1 OP rhs.m1;                                  \
        m2 OP rhs.m2;                                  \
        m3 OP rhs.m3;                                  \
        return *this;                                  \
    }                                                  \
    template <class X>                                 \
    this_type& operator OP(const std::complex<X>& rhs) \
    {                                                  \
        m0 OP rhs.real();                              \
        m1 OP rhs.imag();                              \
        return *this;                                  \
    }

        //----use
        DECLARE_OP_EQUAL(+=)
        DECLARE_OP_EQUAL(-=)
//

#undef DECLARE_OP_EQUAL

        template <class X>
        this_type& operator*=(const std::complex<X>& rhs)
        {
            const T ar = static_cast<T>(rhs.real());
            const T br = static_cast<T>(rhs.imag());

            m0 = (+m0 * ar - m1 * br);
            m1 = (+m0 * br + m1 * ar);
            m2 = (+m2 * ar + m3 * br);
            m3 = (-m2 * br + m3 * ar);

            return (*this);
        }

        template <class X>
        this_type& operator/=(const std::complex<T>& rhs)
        {
            const T ar = static_cast<T>(rhs.real());
            const T br = static_cast<T>(rhs.imag());

            T denominator = T(1) / (ar * ar + br * br);

            m0 = (+m0 * ar + m1 * br) * denominator; //++
            m1 = (-m0 * br + m1 * ar) * denominator; //-+
            m2 = (+m2 * ar - m3 * br) * denominator; //+-
            m3 = (+m2 * br + m3 * ar) * denominator; //++

            return (*this);
        }

#define RHS_W rhs.w()
#define RHS_X rhs.x()
#define RHS_Y rhs.y()
#define RHS_Z rhs.z()

#define LHS_W lhs.w()
#define LHS_X lhs.x()
#define LHS_Y lhs.y()
#define LHS_Z lhs.z()

        template <class X>
        this_type& operator*=(const quaternion<X>& rhs)
        {
            const T tmp[4] = {m0, m1, m2, m3};

            element[W_INDEX] = +tmp[W_INDEX] * RHS_W - tmp[X_INDEX] * RHS_X - tmp[Y_INDEX] * RHS_Y - tmp[Z_INDEX] * RHS_Z;
            element[X_INDEX] = +tmp[W_INDEX] * RHS_X + tmp[X_INDEX] * RHS_W + tmp[Y_INDEX] * RHS_Z - tmp[Z_INDEX] * RHS_Y; //(m0*RHS_X+RHS_W*b)+(m2*RHS_Z-RHS_Y*d);
            element[Y_INDEX] = +tmp[W_INDEX] * RHS_Y - tmp[X_INDEX] * RHS_Z + tmp[Y_INDEX] * RHS_W + tmp[Z_INDEX] * RHS_X; //(m0*RHS_Y+RHS_W*c)+(m3*RHS_X-RHS_Z*b);
            element[Z_INDEX] = +tmp[W_INDEX] * RHS_Z + tmp[X_INDEX] * RHS_Y - tmp[Y_INDEX] * RHS_X + tmp[Z_INDEX] * RHS_W; //(m0*RHS_Z+RHS_W*d)+(m1*RHS_Y-RHS_X*c);

            return *this;
        }

        template <class X>
        this_type& operator/=(const quaternion<X>& rhs)
        {
            const T tmp[4] = {m0, m1, m2, m3};

            T denominator = T(1) / (RHS_W * RHS_W + RHS_X * RHS_X + RHS_Y * RHS_Y + RHS_Z * RHS_Z);

            element[W_INDEX] = (+tmp[W_INDEX] * RHS_W + tmp[X_INDEX] * RHS_X + tmp[Y_INDEX] * RHS_Y + tmp[Z_INDEX] * RHS_Z) * denominator; //((m0*RHS_W+m1*RHS_X+m2*RHS_Y+m3*RHS_Z)/denominator;
            element[X_INDEX] = (-tmp[W_INDEX] * RHS_X + tmp[X_INDEX] * RHS_W - tmp[Y_INDEX] * RHS_Z + tmp[Z_INDEX] * RHS_Y) * denominator; //((RHS_W*b-m0*RHS_X)+(RHS_Y*d-m2*RHS_Z))/denominator;
            element[Y_INDEX] = (-tmp[W_INDEX] * RHS_Y + tmp[X_INDEX] * RHS_Z + tmp[Y_INDEX] * RHS_W - tmp[Z_INDEX] * RHS_X) * denominator; //((RHS_W*c-m0*RHS_Y)+(RHS_Z*b-m3*RHS_X))/denominator;
            element[Z_INDEX] = (-tmp[W_INDEX] * RHS_Z - tmp[X_INDEX] * RHS_Y + tmp[Y_INDEX] * RHS_X + tmp[Z_INDEX] * RHS_W) * denominator; //((RHS_W*d-m0*RHS_Z)+(RHS_X*c-m1*RHS_Y))/denominator;

            return *this;
        }

#undef RHS_W
#undef RHS_X
#undef RHS_Y
#undef RHS_Z
#undef LHS_W
#undef LHS_X
#undef LHS_Y
#undef LHS_Z

        template <class X>
        this_type& operator*=(const X rhs)
        {
            m0 *= rhs;
            m1 *= rhs;
            m2 *= rhs;
            m3 *= rhs;
            return *this;
        }

        template <class X>
        this_type& operator/=(const X rhs)
        {
            m0 /= rhs;
            m1 /= rhs;
            m2 /= rhs;
            m3 /= rhs;
            return *this;
        }

        //-----------------------------------------------
        T& operator[](size_type i)
        {
            return element[i];
        }

        value_type operator[](size_type i) const
        {
            return element[i];
        }

        reference at(size_type i)
        {
            if (c_size <= i)
            {
                throw std::out_of_range(debug());
            }
            return element[i];
        }
        const_reference at(size_type i) const
        {
            if (c_size <= i)
            {
                throw std::out_of_range(debug());
            }
            return element[i];
        }
        //-----------------------------------------------
        T w() const { return element[W_INDEX]; }
        T x() const { return element[X_INDEX]; }
        T y() const { return element[Y_INDEX]; }
        T z() const { return element[Z_INDEX]; }

        //-----------------------------------------------
        //utilities

        T sqr_length() const
        {
            return m0 * m0 + m1 * m1 + m2 * m2 + m3 * m3;
        }
        T length() const
        {
            using namespace std;
            return sqrt(sqr_length());
        }
        T dot(const quaternion<T>& rhs) const
        {
            return m0 * rhs.m0 + m1 * rhs.m1 + m2 * rhs.m2 + m3 * rhs.m3;
        }

        /*
         * euc_
         *
         */
        quaternion<T>& normalize()
        {
            using namespace std;
            T length = sqr_length(); //||q||^2
            //if (length == T()) return *this;

            length = T(1) / sqrt(length); // 1 / ||q||

            m0 *= length;
            m1 *= length;
            m2 *= length;
            m3 *= length;

            return *this;
        }

        /** 
         * @name others
         */
        //@{
        const char* debug() const { return "tempest::quaternion<T>"; }
        const char* order() const { return "wxyz"; }
        //@}
    private:
        union
        {
            struct
            {
                T m0, m1, m2, m3;
            };
            T element[4];
        };
    };

    /**
     * @relates quaternion
     */
    //@{

    template <class T>
    inline quaternion<T> operator+(const quaternion<T>& rhs)
    {
        return rhs;
    }

    template <class T>
    inline quaternion<T> operator-(const quaternion<T>& rhs)
    {
        return quaternion<T>(rhs).negate();
    }
    //-----------------------------------------------
    //binary operators

    template <class T>
    inline quaternion<T> operator+(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        return quaternion<T>(lhs) += rhs;
    }
    template <class T>
    inline quaternion<T> operator-(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        return quaternion<T>(lhs) -= rhs;
    }

    template <class T>
    inline quaternion<T> operator*(T lhs, const quaternion<T>& rhs)
    {
        return quaternion<T>(rhs) *= lhs;
    }

    template <class T>
    inline quaternion<T> operator*(const quaternion<T>& lhs, T rhs)
    {
        return quaternion<T>(lhs) *= rhs;
    }

    template <class T>
    inline quaternion<T> operator/(const quaternion<T>& lhs, T rhs)
    {
        return quaternion<T>(lhs) /= rhs;
    }

#define RHS_W rhs.w()
#define RHS_X rhs.x()
#define RHS_Y rhs.y()
#define RHS_Z rhs.z()

#define LHS_W lhs.w()
#define LHS_X lhs.x()
#define LHS_Y lhs.y()
#define LHS_Z lhs.z()

    template <class T>
    inline quaternion<T> operator*(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        return quaternion<T>::wxyz( //mul4*4 sub*6 add*6 = 28
            +LHS_W * RHS_W - LHS_X * RHS_X - LHS_Y * RHS_Y - LHS_Z * RHS_Z,
            +LHS_W * RHS_X + LHS_X * RHS_W + LHS_Y * RHS_Z - LHS_Z * RHS_Y, //(m0*RHS_X+RHS_W*b)+(m2*RHS_Z-RHS_Y*d);
            +LHS_W * RHS_Y - LHS_X * RHS_Z + LHS_Y * RHS_W + LHS_Z * RHS_X, //(m0*RHS_Y+RHS_W*c)+(LHS_Z*RHS_X-RHS_Z*b);
            +LHS_W * RHS_Z + LHS_X * RHS_Y - LHS_Y * RHS_X + LHS_Z * RHS_W  //(m0*RHS_Z+RHS_W*d)+(m1*RHS_Y-RHS_X*c);
            );
    }

    template <class T>
    inline quaternion<T> operator/(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        const T denominator = T(1) / (RHS_W * RHS_W + RHS_X * RHS_X + RHS_Y * RHS_Y + RHS_Z * RHS_Z);

        return quaternion<T>::wxyz(
            (+LHS_W * RHS_W + LHS_X * RHS_X + LHS_Y * RHS_Y + LHS_Z * RHS_Z) * denominator, //(m0*RHS_W+m1*RHS_X+m2*RHS_Y+m3*RHS_Z)/denominator;
            (-LHS_W * RHS_X + LHS_X * RHS_W - LHS_Y * RHS_Z + LHS_Z * RHS_Y) * denominator, //((RHS_W*b-m0*RHS_X)+(RHS_Y*d-m2*RHS_Z))/denominator;
            (-LHS_W * RHS_Y + LHS_X * RHS_Z + LHS_Y * RHS_W - LHS_Z * RHS_X) * denominator, //((RHS_W*c-m0*RHS_Y)+(RHS_Z*b-m3*RHS_X))/denominator;
            (-LHS_W * RHS_Z - LHS_X * RHS_Y + LHS_Y * RHS_X + LHS_Z * RHS_W) * denominator  //((RHS_W*d-m0*RHS_Z)+(RHS_X*c-m1*RHS_Y))/denominator;
            );
    }

#undef RHS_W
#undef RHS_X
#undef RHS_Y
#undef RHS_Z
#undef LHS_W
#undef LHS_X
#undef LHS_Y
#undef LHS_Z

    //-----------------------------------------------
    // utility functions

    template <class T, std::size_t Sz>
    inline quaternion<T> normalize(const quaternion<T>& rhs)
    {
        return quaternion<T>(rhs).normalize();
    }

    template <class T, std::size_t Sz>
    inline T length(const quaternion<T>& rhs)
    {
        return rhs.length();
    }

    template <class T, std::size_t Sz>
    inline T sqr_length(const quaternion<T>& rhs)
    {
        return rhs.sqr_length();
    }

    template <class T>
    inline T dot(const quaternion<T>& lhs, const quaternion<T>& rhs)
    {
        return lhs.dot(rhs);
    }

    template <class T>
    inline quaternion<T> operator~(const quaternion<T>& rhs)
    {
        typedef quaternion<T> QT;

        return QT::wxyz(rhs.w(), -rhs.x(), -rhs.y(), -rhs.z());
    }

    template <class T>
    inline quaternion<T> conj(const quaternion<T>& rhs)
    { //Conjugate
        typedef quaternion<T> QT;

        T length = rhs.sqr_length(); //lq|*|q|
        //if (l == T()) return rhs;
        length = T(1) / length;
        return QT::wxyz(rhs.w() * length, -rhs.x() * length, -rhs.y() * length, -rhs.z() * length);
    }

    /**
     * q * ~q
     */
    template <class T>
    inline T norm(const quaternion<T>& rhs)
    {
        return rhs.sqr_length();
        //return (rhs * ~rhs)[QT::W_INDEX];
    }

    template <class T>
    inline quaternion<T> lerp(const quaternion<T>& lhs, const quaternion<T>& rhs, T t)
    {
        return ((1 - t) * lhs + t * rhs).normalize();
    }

    template <class T>
    inline quaternion<T> slerp(const quaternion<T>& lhs, const quaternion<T>& rhs, T t)
    {
        using namespace std;
        T theta = acos(dot(lhs, rhs));
        T s = static_cast<T>(T(1) / sin(theta));
        return (sin((1 - t) * theta) * lhs + sin(t * theta) * rhs) * s;
    }

    /*
    template<class X,class Y>
    quaternion(X s, const vector<Y,3> & v):m0(s),m1(v[QT::W_INDEX]),m2(v[QT::X_INDEX]),m3(v[QT::Y_INDEX]){}//scalar , vector
    */
    template <class Self, class T>
    quaternion<T> scalar_vector(T s, const vector_base<Self, T, 3>& bv)
    { //scalar , vector
        const Self& v = static_cast<const Self&>(bv);
        return quaternion<T>::wxyz(s, v[0], v[1], v[2]);
    }

    template <class Self, class T>
    quaternion<T> pivot_angle(const vector_base<Self, T, 3>& bv, T theta)
    { //pivot , theta
        using namespace std;
        const Self& v = static_cast<const Self&>(bv);

        theta /= T(2);
        T s = sin(theta);
		T c = cos(theta);
        return quaternion<T>::wxyz(c, s * v[0], s * v[1], s * v[2]);
    }

    template <typename T, typename _CharT, class _Traits>
    std::basic_ostream<_CharT, _Traits>& operator<<(std::basic_ostream<_CharT, _Traits>& os, const quaternion<T>& rhs)
    {
        std::basic_ostringstream<_CharT, _Traits> s;
        s.flags(os.flags());
        s.imbue(os.getloc());
        s.precision(os.precision());
#if 1
        s << "(" << rhs[0] << "," << rhs[1] << "," << rhs[2] << "," << rhs[3] << ")";
#else
        s << "(" << rhs.w() << "," << rhs.x() << "," << rhs.y() << "," << rhs.z() << ")";
#endif
        return os << s.str();
    }

    //@}

    //-----------------------------------------------
    //template<>    class quaternion<float>;
    //template<double>    class quaternion<double>;
    //template<long double>    class quaternion<long double>;

} //

#endif
