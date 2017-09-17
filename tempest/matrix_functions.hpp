#pragma once
#ifndef TM_MATRIX_FUNCTIONS_HPP
#define TM_MATRIX_FUNCTIONS_HPP
/**
 * @file matrix_functions.hpp
 * @brief Defines matrix functions.
 *
 * $RCSfile: matrix_functions.hpp,v $
 * $Date: 2005/02/26 14:39:01 $
 * $Author: Toru Matsuoka $
 */

//#include<iosfwd>
//#include<ostream>
#include <sstream>

#include <algorithm> //fill_n

namespace tempest
{

    //-----------------------------------------------
    //Not a member!
    //-----------------------------------------------

    //-----------------------------------------------
    //operators

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline matrix<T, RowSz, ColumnSz> operator+(const matrix<T, RowSz, ColumnSz>& rhs)
    {
        return rhs;
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline matrix<T, RowSz, ColumnSz> operator-(const matrix<T, RowSz, ColumnSz>& rhs)
    {
        return matrix<T, RowSz, ColumnSz>(rhs).negate();
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline matrix<T, RowSz, ColumnSz> operator+(const matrix<T, RowSz, ColumnSz>& lhs, const matrix<T, RowSz, ColumnSz>& rhs)
    {
        return matrix<T, RowSz, ColumnSz>(lhs) += rhs;
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline matrix<T, RowSz, ColumnSz> operator-(const matrix<T, RowSz, ColumnSz>& lhs, const matrix<T, RowSz, ColumnSz>& rhs)
    {
        return matrix<T, RowSz, ColumnSz>(lhs) -= rhs;
    }

    template <class T, std::size_t ASz, std::size_t BSz, std::size_t CSz> //This functon is NOT recommended.
    matrix<T, ASz, CSz> operator*(const matrix<T, ASz, BSz>& lhs, const matrix<T, BSz, CSz>& rhs)
    {
        matrix<T, ASz, CSz> tmp;
        std::fill_n(tmp.begin(), ASz * CSz, T());

        for (std::size_t i = 0; i < ASz; ++i)
        {
            for (std::size_t k = 0; k < BSz; k++)
            {
                for (std::size_t j = 0; j < CSz; ++j)
                {
                    tmp[i][j] += lhs[i][k] * rhs[k][j];
                }
            }
        }
        return tmp;
    }

    template <class T, std::size_t ASz, std::size_t BSz, std::size_t CSz>
    void multiply(matrix<T, ASz, CSz>* out, const matrix<T, ASz, BSz>& lhs, const matrix<T, BSz, CSz>& rhs)
    {
        matrix<T, ASz, CSz> tmp;
        std::fill_n<T>((*out).begin(), ASz * CSz, T());

        for (std::size_t i = 0; i < ASz; ++i)
        {
            for (std::size_t k = 0; k < BSz; k++)
            {
                for (std::size_t j = 0; j < CSz; ++j)
                {
                    tmp[i][j] += lhs[i][k] * rhs[k][j];
                }
            }
        }

        *out = tmp;
    }

    //-----------------------------------------------
    //specific scalar

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline matrix<T, RowSz, ColumnSz> operator*(T lhs, const matrix<T, RowSz, ColumnSz>& rhs)
    {
        return matrix<T, RowSz, ColumnSz>(rhs) *= lhs;
    }
    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline matrix<T, RowSz, ColumnSz> operator*(const matrix<T, RowSz, ColumnSz>& lhs, T rhs)
    {
        return matrix<T, RowSz, ColumnSz>(lhs) *= rhs;
    }
    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline matrix<T, RowSz, ColumnSz> operator/(T lhs, const matrix<T, RowSz, ColumnSz>& rhs)
    {
        return matrix<T, RowSz, ColumnSz>(rhs) /= lhs;
    }
    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline matrix<T, RowSz, ColumnSz> operator/(const matrix<T, RowSz, ColumnSz>& lhs, T rhs)
    {
        return matrix<T, RowSz, ColumnSz>(lhs) /= rhs;
    }

    //-----------------------------------------------
    // utilities
    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    matrix<T, ColumnSz, RowSz> transpose(const matrix<T, RowSz, ColumnSz>& rhs)
    {

        matrix<T, ColumnSz, RowSz> tmp;
        for (std::size_t i = 0; i < ColumnSz; i++)
        {
            for (std::size_t j = 0; j < RowSz; j++)
            {
                tmp[i][j] = rhs[j][i];
            }
        }
        return tmp;
    }

    //--------------------------------------------------
    //compare

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    bool operator==(const matrix<T, RowSz, ColumnSz>& lhs, const matrix<T, RowSz, ColumnSz>& rhs)
    {
        typename matrix<T, RowSz, ColumnSz>::const_iterator l_i = lhs.begin();
        typename matrix<T, RowSz, ColumnSz>::const_iterator r_i = rhs.begin();

        typename matrix<T, RowSz, ColumnSz>::const_iterator end = lhs.end();
        while (l_i != end)
        {
            if (*l_i != *r_i)
                return false;
            ++l_i;
            ++r_i;
        }
        return true;
    }

    template <class T, std::size_t RowSz, std::size_t ColumnSz>
    inline bool operator!=(const matrix<T, RowSz, ColumnSz>& lhs, const matrix<T, RowSz, ColumnSz>& rhs)
    {
        return !(lhs == rhs);
    }

/*
template<class T,std::size_t Sz>
inline bool operator< (const matrix<T,Sz> &lhs,const matrix<T,Sz> &rhs){
    return lhs.sqr_length() < rhs.sqr_length();
}

template<class T,std::size_t Sz>
inline bool operator> (const matrix<T,Sz> &lhs,const matrix<T,Sz> &rhs){
    return rhs < lhs;
}

template<class T,std::size_t Sz>
inline bool operator>= (const matrix<T,Sz> &lhs,const matrix<T,Sz> &rhs){
    return !(lhs < rhs);
}

template<class T,std::size_t Sz>
inline bool operator<= (const matrix<T,Sz> &lhs,const matrix<T,Sz> &rhs){
    return !(rhs < lhs);
}
*/

    template <class T, std::size_t Sz>
    bool lu_decomposite(matrix<T, Sz, Sz>& Out, std::size_t indices[Sz], const matrix<T, Sz, Sz>& A)
    {
        int rows = Sz;
        int cols = Sz;

        Out = A;

        for (int k = 0; k < rows; k++)
        {
            indices[k] = k;
        }

        for (int k = 0; k < rows - 1; k++)
        {
            {
                T mkk = std::fabs(Out[k][k]);
                int imax = k;
                for (int i = k + 1; i < rows; i++)
                {
                    T mik = std::fabs(Out[i][k]);
                    if (mik > mkk)
                    {
                        mkk = mik;
                        imax = i;
                    }
                }
                if (mkk == T(0))
                    return false;
                if (imax != k)
                {
                    std::swap(indices[k], indices[imax]);
                }
            }

            {
                int kk = indices[k];
                T mkk = T(1) / Out[kk][k];
                for (int i = k + 1; i < rows; i++)
                {
                    int ii = indices[i];
                    Out[ii][k] *= mkk;
                }
            }
            {
                int kk = indices[k];
                for (int i = k + 1; i < rows; i++)
                {
                    for (int j = k + 1; j < cols; j++)
                    {
                        int ii = indices[i];
                        Out[ii][j] -= Out[ii][k] * Out[kk][j];
                    }
                }
            }
        }

        return true;
    }
    template <class T, std::size_t Sz> //lu
    bool invert_lu(matrix<T, Sz, Sz>& Out, const matrix<T, Sz, Sz>& A)
    {
        int rows = Sz;
        int cols = Sz;

        matrix<T, Sz, Sz> B(A);
        std::size_t indices[Sz];
        if (!lu_decomposite(B, indices, A))
            return false;
        matrix<T, Sz, Sz>& C = Out;
        C = B;

        for (int k = 0; k < rows; k++)
        {
            // lower
            for (int i = 0; i < rows; i++)
            {
                T t = T(i == k ? 1 : 0);
                for (int j = 0; j < i; j++)
                {
                    t -= B[indices[i]][j] * C[j][k];
                }
                C[i][k] = t;
            }

            // upper
            for (int i = rows - 1; i >= 0; i--)
            {
                T t = C[indices[i]][k];
                for (int j = i + 1; j < cols; j++)
                {
                    t -= B[indices[i]][j] * C[j][k];
                }
                C[i][k] = t / B[indices[i]][i];
            }
        }

        return true;
    }

    template <class T, std::size_t Sz> //lu
    bool invert_lu(matrix<T, Sz, Sz>* Out, const matrix<T, Sz, Sz>& A)
    {
        if (Out == 0)
            return false;
        return invert_lu(*Out, A);
    }

    template <class T, std::size_t Sz> //gauss joldan
    bool invert_gj(matrix<T, Sz, Sz>* out, const matrix<T, Sz, Sz>& rhs)
    {
        int i, j, k;
        T t, u, det;

        if (out == 0)
            return false;

        matrix<T, Sz, Sz>& a = static_cast<matrix<T, Sz, Sz>&>(*out);

        a = rhs;

        det = T(1);
        for (k = 0; k < Sz; k++)
        {
            t = a[k][k];
            det *= t;
            for (i = 0; i < Sz; i++)
                a[k][i] /= t;
            a[k][k] = T(1) / t;
            for (j = 0; j < Sz; j++)
            {
                if (j != k)
                {
                    u = a[j][k];
                    for (i = 0; i < Sz; i++)
                        if (i != k)
                            a[j][i] -= a[k][i] * u;
                        else
                            a[j][i] = -u / t;
                }
            }
        }

        if (det != T(0))
            return true;
        else
            return false;
        //return det;
    }

    //if true definition matrix invert function using by LU separate method.

    //#define __TM_INVERT_USE_LU__

    template <class T, std::size_t Sz> //gauss joldan
    inline bool invert(matrix<T, Sz, Sz>* out, const matrix<T, Sz, Sz>& rhs)
    {
#ifdef __TM_INVERT_USE_LU__
        return invert_lu(out, rhs);
#else
        return invert_gj(out, rhs);
#endif
    }

    template <class T, std::size_t Sz> //gauss joldan
    inline matrix<T, Sz, Sz> invert(const matrix<T, Sz, Sz>& rhs)
    {
        matrix<T, Sz, Sz> tmp;
        bool bRet = invert(&tmp, rhs);
        assert(bRet);
        return tmp;
    }

    //

    //-----------------------------------------------
    //output

    /**
     * ostream <<
     */
    template <typename T, std::size_t RowSz, std::size_t ColumnSz, typename _CharT, class _Traits>
    std::basic_ostream<_CharT, _Traits>& operator<<(std::basic_ostream<_CharT, _Traits>& os, const matrix<T, RowSz, ColumnSz>& rhs)
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

} //end of namespace

#endif
