#pragma once
#ifndef TM_MATRIX_GENERATOR_HPP
#define TM_MATRIX_GENERATOR_HPP
/**
 *  @file matrix_generator.hpp
 *  @brief Defines class template tempest::matrix_generator<T>,which is helper for generating matrix.
 *  @author Toru Matsuoka
 *  $Id: matrix_generator.hpp,v 1.5 2005/03/06 14:27:55      Exp $
 */
#include <cmath>
#include <cstddef>
#include <complex>

namespace tempest
{

    namespace ns_mg
    {
        template <class T>
        inline T sin(T x)
        {
            return std::sin(x);
        }
        template <class T>
        inline T cos(T x)
        {
            return std::cos(x);
        }
    }

/**
 *  @struct matrix_generator
 *  @brief  Helper class template for matrix generation.
 *
 *  @code
 *  typedef tempest::matrix<double,2> mat2;                         //matrix<double,2>
 *  typedef tempest::matrix_generator< mat2 > mgen2;  //generator for matrix<double,2>
 *  mat2 m(mgen2::scaling(1.5,1.5));
 *  @endcode
 *
 */

/**
 *  @struct matrix_generator< matrix<T, Sz> >
 *  @brief  Helper class template for matrix generation.
 *
 *  @code
 *  typedef tempest::matrix<double,10> mat10;                         //matrix<double,10>
 *  typedef tempest::matrix_generator< mat10 > mgen10;  //generator for matrix<double,10>
 *  mat10 m(mgen10::zero());                                //0
 *  @endcode
 */

/**
 *  @struct matrix_generator< matrix<T,2,2> >
 *  @brief  Helper class template for matrix generation.
 *
 *  @code
 *  typedef tempest::matrix<double,2> mat2;                         //matrix<double,2>
 *  typedef tempest::matrix_generator< mat2 > mgen2;  //generator for matrix<double,2>
 *  mat2 m(mgen2::scaling(1.5,1.5));
 *  @endcode
 */

/**
 *  @struct matrix_generator< matrix<T,3,3> >
 *  @brief  Helper class template for matrix generation.
 *
 *  @code
 *  typedef tempest::matrix<double,3> mat3;                         //matrix<double,3>
 *  typedef tempest::matrix_generator< mat3 > mgen3;  //generator for matrix<double,3>
 *  mat3 m(mgen3::scaling(1.5,1.5,1.5));
 *  @endcode
 *
 */

/**
 *  @struct matrix_generator< matrix<T,4,4> >
 *  @brief  Helper class template for matrix generation.
 *
 *  @code
 *  typedef tempest::matrix<double,4> mat4;                         //matrix<double,4>
 *  typedef tempest::matrix_generator< mat4 > mgen4;  //generator for matrix<double,4>
 *  mat4 m(mgen4::scaling(1.5,1.5,1.5));
 *  @endcode
 */

    template <class T>
    struct matrix_generator
    {
    };

    template <class T, std::size_t Sz>
    struct matrix_generator<matrix<T, Sz, Sz> >
    {
        static matrix<T, Sz, Sz> zero()
        {
            matrix<T, Sz> tmp;
            for (std::size_t i = 0; i < Sz; i++)
            {
                for (std::size_t j = 0; j < Sz; j++)
                {
                    tmp[i][j] = T();
                }
            }

            return tmp;
        }

        static matrix<T, Sz, Sz> identity()
        {
            matrix<T, Sz, Sz> tmp;
            for (std::size_t i = 0; i < Sz; i++)
            {
                for (std::size_t j = 0; j < Sz; j++)
                {
                    tmp[i][j] = T(0);
                }
            }

            for (std::size_t i = 0; i < Sz; i++)
            {
                tmp[i][i] = T(1);
            }

            return tmp;
        }
    };

    //-------------------------------------------------------------------
    //For "matrix2"
    template <class T>
    struct matrix_generator<matrix<T, 2, 2> >
    {
        static inline matrix<T, 2, 2> zero()
        {
            return matrix<T, 2, 2>(0, 0, 0, 0);
        }

        static inline matrix<T, 2, 2> identity()
        {
            return matrix<T, 2, 2>(1, 0, 0, 1);
        }

        static inline matrix<T, 2, 2> scaling(const T _x, const T _y)
        {
            return matrix<T, 2, 2>(_x, 0, 0, _y);
        }

        static inline matrix<T, 2, 2> rotation_z(const T s, const T c)
        {
            return matrix<T, 2, 2>(c, -s, s, c);
        }

        static inline matrix<T, 2, 2> rotation_z(const T radian)
        {
            T _sin = ns_mg::sin<T>(radian);
            T _cos = ns_mg::cos<T>(radian);

            return matrix<T, 2, 2>(_cos, -_sin, _sin, _cos);
        }

        static inline matrix<T, 2, 2> convert(const std::complex<T>& rhs)
        {
            T _r = rhs.real();
            T _i = rhs.imag();

            return matrix<T, 2, 2>(
                _r, -_i,
                _i, _r);
        }
    };

    //-------------------------------------------------------------------
    //For matrix3
    template <class T>
    struct matrix_generator<matrix<T, 3, 3> >
    {

        static inline matrix<T, 3, 3> zero()
        {
            return matrix<T, 3, 3>(
                0, 0, 0,
                0, 0, 0,
                0, 0, 0);
        }

        static inline matrix<T, 3, 3> identity()
        {
            return matrix<T, 3, 3>(
                1, 0, 0,
                0, 1, 0,
                0, 0, 1);
        }

        static inline matrix<T, 3, 3> translation(const T tx, const T ty)
        {
            return matrix<T, 3, 3>(
                1, 0, tx,
                0, 1, ty,
                0, 0, 1);
        }

        static inline matrix<T, 3, 3> scaling(const T sx, const T sy, const T sz)
        {
            return matrix<T, 3, 3>(
                sx, 0, 0,
                0, sy, 0,
                0, 0, sz);
        }

        static inline matrix<T, 3, 3> rotation_z(const T s, const T c)
        {
            return matrix<T, 3, 3>(
                c, -s, 0,
                s, c, 0,
                0, 0, 1);
        }

        static inline matrix<T, 3, 3> rotation_z(const T radian)
        {
            T _sin = ns_mg::sin<T>(radian);
            T _cos = ns_mg::cos<T>(radian);
            return rotation_z(_sin, _cos);
        }
        //-----------------------------------------------
        static inline matrix<T, 3, 3> rotation_x(const T s, const T c)
        {
            return matrix<T, 3, 3>(
                1, 0, 0,
                0, c, -s,
                0, s, c);
        }

        static inline matrix<T, 3, 3> rotation_x(const T radian)
        {
            T _sin = static_cast<T>(ns_mg::sin<T>(radian));
            T _cos = static_cast<T>(ns_mg::cos<T>(radian));
            return rotation_x(_sin, _cos);
        }

        //-----------------------------------------------
        static inline matrix<T, 3, 3> rotation_y(const T s, const T c)
        {
            return matrix<T, 3, 3>(
                c, 0, s,
                0, 1, 0,
                -s, 0, c);
        }

        static inline matrix<T, 3, 3> rotation_y(const T radian)
        {
            T _sin = static_cast<T>(ns_mg::sin<T>(radian));
            T _cos = static_cast<T>(ns_mg::cos<T>(radian));
            return rotation_y(_sin, _cos);
        }

        static matrix<T, 3, 3> convert(const quaternion<T>& rhs)
        {
            //4C2 = 4*5/2*1 = 10
            //tww = 1;
            //10 - 1 = 9
            //T tww = 2*rhs[0]*rhs[0];//0
            T w = rhs.w();
            T x = rhs.x();
            T y = rhs.y();
            T z = rhs.z();

            T twx = 2 * w * x; //1
            T twy = 2 * w * y; //2
            T twz = 2 * w * z; //3

            T txx = 2 * x * x; //4
            T txy = 2 * x * y; //5
            T txz = 2 * x * z; //6

            T tyy = 2 * y * y; //7
            T tyz = 2 * z * y; //8

            T tzz = 2 * z * z; //9

            return matrix<T, 3, 3>(
                1 - tyy - tzz, txy - twz, txz + twy,
                txy + twz, 1 - txx - tzz, tyz - twx,
                txz - twy, tyz + twx, 1 - txx - tyy);
        }
    };

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    //For matrix4

    template <class T>
    struct matrix_generator<matrix<T, 4, 4> >
    {
        //typedef typename matrix_generator< matrix<T,4,4> > this_type;

        static inline matrix<T, 4, 4> zero()
        {
            return matrix<T, 4, 4>(
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0);
        }

        static inline matrix<T, 4, 4> identity()
        {
            return matrix<T, 4, 4>(
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1);
        }

        static inline matrix<T, 4, 4> translation(const T tx, const T ty, const T tz)
        {
            return matrix<T, 4, 4>(
                1, 0, 0, tx,
                0, 1, 0, ty,
                0, 0, 1, tz,
                0, 0, 0, 1);
        }

        static inline matrix<T, 4, 4> scaling(const T sx, const T sy, const T sz)
        {
            return matrix<T, 4, 4>(
                sx, 0, 0, 0,
                0, sy, 0, 0,
                0, 0, sz, 0,
                0, 0, 0, 1);
        }

        static inline matrix<T, 4, 4> rotation_z(const T s, const T c)
        {
            return matrix<T, 4, 4>(
                c, -s, 0, 0,
                s, c, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1);
        }

        static inline matrix<T, 4, 4> rotation_z(const T radian)
        {
            T _sin = ns_mg::sin<T>(radian);
            T _cos = ns_mg::cos<T>(radian);
            return rotation_z(_sin, _cos);
        }
        //-----------------------------------------------
        static inline matrix<T, 4, 4> rotation_x(const T s, const T c)
        {
            return matrix<T, 4, 4>(
                1, 0, 0, 0,
                0, c, -s, 0,
                0, s, c, 0,
                0, 0, 0, 1);
        }

        static inline matrix<T, 4, 4> rotation_x(const T radian)
        {
            T _sin = ns_mg::sin<T>(radian);
            T _cos = ns_mg::cos<T>(radian);
            return rotation_x(_sin, _cos);
        }

        //-----------------------------------------------
        static inline matrix<T, 4, 4> rotation_y(const T s, const T c)
        {
            return matrix<T, 4, 4>(
                c, 0, s, 0,
                0, 1, 0, 0,
                -s, 0, c, 0,
                0, 0, 0, 1);
        }

        static inline matrix<T, 4, 4> rotation_y(const T radian)
        {
            T _sin = ns_mg::sin<T>(radian);
            T _cos = ns_mg::cos<T>(radian);
            return rotation_y(_sin, _cos);
        }

        static inline matrix<T, 4, 4> convert(const quaternion<T>& rhs)
        {
            //typedef quaternion<T> QT;
            //4C2 = 4*5/2*1 = 10
            //tww = 1;
            //10 - 1 = 9
            //T tww = 2*rhs[0]*rhs[0];//0
            T w = rhs.w();
            T x = rhs.x();
            T y = rhs.y();
            T z = rhs.z();

            T twx = 2 * w * x; //1
            T twy = 2 * w * y; //2
            T twz = 2 * w * z; //3

            T txx = 2 * x * x; //4
            T txy = 2 * x * y; //5
            T txz = 2 * x * z; //6

            T tyy = 2 * y * y; //7
            T tyz = 2 * z * y; //8

            T tzz = 2 * z * z; //9

            return matrix<T, 4, 4>(
                1 - tyy - tzz, txy - twz, txz + twy, 0,
                txy + twz, 1 - txx - tzz, tyz - twx, 0,
                txz - twy, tyz + twx, 1 - txx - tyy, 0,
                0, 0, 0, 1);
        }
    };

} //end of namespace tempest;
#endif
