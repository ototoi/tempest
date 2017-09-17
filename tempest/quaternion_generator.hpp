#pragma once
#ifndef TM_QUATERNION_GENERATOR_HPP
#define TM_QUATERNION_GENERATOR_HPP

/** 
 *  @file quaternion_generator.hpp
 *  @brief Defines tempest::quaternion_generator< tempest::quaternion<T> >.
 *  @author Toru Matsuoka
 */

#include <cstddef>
#include <cmath>
#include "vector.hpp"
#include "vector3.hpp"
#include "matrix.hpp"
#include "matrix3.hpp"
#include "matrix4.hpp"
#include "quaternion.hpp"

namespace tempest
{
    template <class X>
    struct quaternion_generator
    {
    };

    template <class T>
    struct quaternion_generator< quaternion<T> >
    {
        static inline quaternion<T> convert(const matrix<T, 3, 3>& m)
        {
            return from_matrix(m);
        }

        static inline quaternion<T> convert(const matrix<T, 4, 4>& m)
        {
            return from_matrix(m);
        }
        //---------------------------------------------------------------
#define m00 m[0][0]
#define m10 m[1][0]
#define m20 m[2][0]
#define m30 m[3][0]

#define m01 m[0][1]
#define m11 m[1][1]
#define m21 m[2][1]
#define m31 m[3][1]

#define m02 m[0][2]
#define m12 m[1][2]
#define m22 m[2][2]
#define m32 m[3][2]

#define m03 m[0][3]
#define m13 m[1][3]
#define m23 m[2][3]
#define m33 m[3][3]
        static inline quaternion<T> from_matrix(const matrix<T, 3, 3>& m)
        {
            T tr = m00 + m11 + m22;
    
            T qw, qx, qy, qz;
            if (tr > T(0))
            {
                T S = sqrt(tr + T(1.0)) * T(2.0); // S=4*qw 
                qw = T(0.25) * S;
                qx = (m21 - m12) / S;
                qy = (m02 - m20) / S;
                qz = (m10 - m01) / S;
            }
            else if ((m00 > m11) & (m00 > m22))
            {
                T S = sqrt(T(1.0) + m00 - m11 - m22) * T(2.0); // S=4*qx 
                qw = (m21 - m12) / S;
                qx = T(0.25) * S;
                qy = (m01 + m10) / S;
                qz = (m02 + m20) / S;
            }
            else if (m11 > m22)
            {
                T S = sqrt(T(1.0) + m11 - m00 - m22) * T(2.0); // S=4*qy
                qw = (m02 - m20) / S;
                qx = (m01 + m10) / S;
                qy = T(0.25) * S;
                qz = (m12 + m21) / S;
            }
            else
            {
                T S = sqrt(T(1.0) + m22 - m00 - m11) * T(2.0); // S=4*qz
                qw = (m10 - m01) / S;
                qx = (m02 + m20) / S;
                qy = (m12 + m21) / S;
                qz = T(0.25) * S;
            }
            return quaternion<T>::wxyz(qw, qx, qy, qz);
        }
#undef m00
#undef m10
#undef m20
#undef m30

#undef m01
#undef m11
#undef m21
#undef m31

#undef m02
#undef m12
#undef m22
#undef m32

#undef m03
#undef m13
#undef m23
#undef m33

        static inline quaternion<T> from_matrix(const matrix<T, 4, 4>& m)
        {
            matrix<T, 3, 3> m3;
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    m3[i][j] = m[i][j];
                }
            }
            return from_matrix(m3);
        }

        //---------------------------------------------------------------

        static inline quaternion<T> pivot_angle(const vector<T, 3>& v, T r)
        {
            return tempest::pivot_angle(v, r);
        }

        static inline quaternion<T> from_pivot_angle(const vector<T, 3>& v, T r)
        {
            return pivot_angle(v, r);
        }

        //---------------------------------------------------------------

        static inline quaternion<T> between(const vector<T, 3>& ax, const vector<T, 3>& bx)
        {
            vector<T, 3> an = normalize(ax);
            vector<T, 3> bn = normalize(bx);
            vector<T, 3> pivot = cross(an, bn);
            T c = dot(an, bn);
            T l = pivot.length();
            if(l > T(0))
            {
                T angle = asin(l);//0 ~ +PI/2
                if (c < T(0))
                {
                    angle = pi() - angle;
                }
                pivot *= (T(1) / l);
                return pivot_angle(pivot, angle);
            }
            else
            {
                return quaternion<T>::wxyz(1, 0, 0, 0);
            }
        }

        static inline quaternion<T> from_two_axis(const vector<T, 3>& ax, const vector<T, 3>& bx)
        {
            return between(ax, bx);
        }

        //---------------------------------------------------------------

        static inline vector<T, 3> quaternion_to_exponentialmap(const quaternion<T>& q)
        {
            T theta = acos(q.w());
            vector<T, 3> v(q.x(), q.y(), q.z());
            if(v.length() > T(0))
            {
                v.normalize();
            }
            return v * theta;
        }

        static inline quaternion<T> exponentialmap_to_quaternion(const vector<T, 3>& v)
        {
            T theta = v.length();
            vector<T, 3> v2 = v;
            if(v2.length() > T(0))
            {
                v2.normalize();
            }
            T c = cos(theta);
            T s = sin(theta);
            return quaternion<T>::wxyz(c, s*v2[0], s*v2[1], s*v2[2]);
        }

        static inline quaternion<T> from_exponentialmap(const vector<T, 3>& v)
        {
            return exponentialmap_to_quaternion(v);
        }

        //---------------------------------------------------------------

        static inline quaternion<T> from_euler_angles(T rx, T ry, T rz)
        {
            return from_euler_angles_xyz(rx, ry, rz);
        }

        //---------------------------------
        static inline quaternion<T> from_euler_angles_xyz(T rx, T ry, T rz)
        {
            return ax(rx) * ay(ry) * az(rz);
        }

        static inline quaternion<T> from_euler_angles_xzy(T rx, T ry, T rz)
        {
            return ax(rx) * az(rz) * ay(ry);
        }

        static inline quaternion<T> from_euler_angles_yzx(T rx, T ry, T rz)
        {
            return ay(ry) * az(rz) * ax(rx);
        }

        static inline quaternion<T> from_euler_angles_yxz(T rx, T ry, T rz)
        {
            return ay(ry) * ax(rx) * az(rz);
        }

        static inline quaternion<T> from_euler_angles_zxy(T rx, T ry, T rz)
        {
            return az(rz) * ax(rx) * ay(ry);
        }

        static inline quaternion<T> from_euler_angles_zyx(T rx, T ry, T rz)
        {
            return az(rz) * ay(ry) * ax(rx);
        }

        //---------------------------------
        static inline quaternion<T> rotation_x(T r)
        {
            return ax(r);
        }

        static inline quaternion<T> rotation_y(T r)
        {
            return ay(r);
        }

        static inline quaternion<T> rotation_z(T r)
        {
            return az(r);
        }

        static inline quaternion<T> ax(T r)
        {
            typedef vector<T, 3> V;
            return pivot_angle(V(1, 0, 0), r);
        }

        static inline quaternion<T> ay(T r)
        {
            typedef vector<T, 3> V;
            return pivot_angle(V(0, 1, 0), r);
        }

        static inline quaternion<T> az(T r)
        {
            typedef vector<T, 3> V;
            return pivot_angle(V(0, 0, 1), r);
        }

        //---------------------------------
        static inline T pi()
        {
            return T(3.14159265358979323846);
        }

        static inline T sin(T x)
        {
            return std::sin(x);
        }

        static inline T cos(T x)
        {
            return std::cos(x);
        }

        static inline T asin(T x)
        {
            return std::asin(x);
        }

        static inline T acos(T x)
        {
            return std::acos(x);
        }

        //---------------------------------
        static inline T angle(const quaternion<T>& q)
        {
            //T w = q.w();
            //T x = q.x();
            //T y = q.y();
            //T z = q.z();
            //T ac = w;
            //T as = sqrt(x*x + y*y + z*z);
            //T s = acos(q.w());
            return T(2) * acos(q.w());
        }

        static inline vector<T, 3> pivot(const quaternion<T>& q)
        {
            using namespace std;
            T x = q.x();
            T y = q.y();
            T z = q.z();
            T as = sqrt(x*x + y*y + z*z);
            if(as > T(0))
            {
                T l = T(1) / as;
                return vector<T, 3>(x * l, y * l, z * l);
            }
            else
            {
                return vector<T, 3>(x, y, z);
            }
        }
    };
}

#endif
