#ifndef TM_VECTOR_FUNCTIONS_HPP
#define TM_VECTOR_FUNCTIONS_HPP

/**
 * @file vector_functions.hpp
 * @brief Defines vecor functions.
 *
 * $RCSfile: vector_functions.hpp,v $
 * $Date: 2005/02/26 14:39:01 $
 * $Author: Toru Matsuoka $
 */

//#include<iofwd>
#include <algorithm> //equal
#include <numeric>
//#include<ostream>
#include <sstream>

#include "vector_base.hpp"

namespace tempest
{

    //-----------------------------------------------
    //Not a member!
    //-----------------------------------------------

    //-----------------------------------------------
    //unary operator
    /**
     *  @name unary operator
     *  @relates vector
     */
    //@{

    template <class T, std::size_t Sz>
    inline const vector<T, Sz> operator+(const vector<T, Sz>& rhs)
    {
        return rhs;
    }

    template <class T, std::size_t Sz>
    inline vector<T, Sz> operator-(const vector<T, Sz>& rhs)
    {
        return vector<T, Sz>(rhs).negate();
    }
//@}

//-----------------------------------------------
//binary operator
/**
 * @name binary operator
 * @relates vector
 */
//@{

#define DECLARE_OPERATOR(OP)                                        \
    \
template<class T, std::size_t Sz> \
inline vector<T, Sz>                                                \
    operator OP(const vector<T, Sz>& lhs, const vector<T, Sz>& rhs) \
    {                                                               \
        return vector<T, Sz>(lhs) OP## = rhs;                       \
    \
}
    DECLARE_OPERATOR(+)
    DECLARE_OPERATOR(-)
    DECLARE_OPERATOR(*)
    DECLARE_OPERATOR(/ )

#undef DECLARE_OPERATOR

    //@}

    //-----------------------------------------------
    //specific scalar
    /**
     *    @name specific scalar
     *    @relates vector
     */
    //@{

    /**
     *    @param lhs an any type.
     *    @param rhs a vector.
     *    @return lhs * rhs
     */
    template <class T, std::size_t Sz>
    inline vector<T, Sz> operator*(const T lhs, const vector<T, Sz>& rhs)
    {
        return vector<T, Sz>(rhs) *= T(lhs);
    }
    /**
     *    @param lhs a vector.
     *     @param rhs an any type.
     *    @return lhs * rhs
     */
    template <class T, std::size_t Sz>
    inline vector<T, Sz> operator*(const vector<T, Sz>& lhs, const T rhs)
    {
        return vector<T, Sz>(lhs) *= T(rhs);
    }

    /**
     *    @param lhs a vector.
     *    @param rhs an any type.
     *    @return lhs / rhs
     */
    template <class T, std::size_t Sz>
    inline vector<T, Sz> operator/(const vector<T, Sz>& lhs, const T rhs)
    {
        return vector<T, Sz>(lhs) /= T(rhs);
    }
    //@}

    //-----------------------------------------------
    // utility functions
    /**
     *    @name utility
     *    @relates vector
     */
    //@{

    /**
     *    @param rhs a vector.
     *    @return || rhs ||
     */
    template <class T, std::size_t Sz>
    inline T length(const vector<T, Sz>& rhs)
    {
        return rhs.length();
    }
    /**
     *    @param rhs a vector.
     *    @return || rhs ||^2
     */
    template <class T, std::size_t Sz>
    inline T sqr_length(const vector<T, Sz>& rhs)
    {
        return rhs.sqr_length();
    }

    /**
     *    @param rhs a vector.
     *    @return sigma (rhs)
     */
    template <class T, std::size_t Sz>
    inline T sum(const vector<T, Sz>& rhs)
    {
        return rhs.sum();
    }

    /**
     *    @param rhs a vector.
     *    @return rhs/|| rhs ||
     */
    template <class T, std::size_t Sz>
    inline vector<T, Sz> normalize(const vector<T, Sz>& rhs)
    {
        return vector<T, Sz>(rhs).normalize();
    }

    /*
     *    @param lhs a vector.
     *    @param rhs a vector.
     *    @return lhs &dot; rhs
     */
    template <class T, std::size_t Sz>
    inline T dot(const vector<T, Sz>& lhs, const vector<T, Sz>& rhs)
    {
        /*
        T temp();
        for(std::size_t i = 0;i < Sz;i++){temp += (lhs[i]*rhs[i]);}
        return temp;
        */
        return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), T());
    }
    //@}

    //--------------------------------------------------
    //compare
    /**
     *    @name comparer
     *    @relates vector
     */
    //@{

    /**
     *    @param lhs a vector.
     *    @param rhs a vector.
     *    @return lhs == rhs
     */
    template <class T, std::size_t Sz>
    inline bool operator==(const vector<T, Sz>& lhs, const vector<T, Sz>& rhs)
    {
        return std::equal(&(lhs[0]), &(lhs[0]) + Sz, &(rhs[0]));
        //for(std::size_t i = 0;i<Sz;++i){if(lhs[i]!=rhs[i])return false;}
        //return true;
        //return !strcmp(&(lhs[0]),&(rhs[0])); ///< with pod type version.
    }

    /**
     *    @param lhs a vector.
     *    @param rhs a vector.
     *    @return lhs != rhs
     */
    template <class T, std::size_t Sz>
    inline bool operator!=(const vector<T, Sz>& lhs, const vector<T, Sz>& rhs)
    {
        return !(lhs == rhs);
    }

#if 0
/**
 *    @param lhs a vector.
 *    @param rhs a vector.
 *    @return lhs < rhs
 */
template<class T,std::size_t Sz>
inline bool operator< (const vector<T,Sz> &lhs, const vector<T,Sz> &rhs){
    return lhs.sqr_length() < rhs.sqr_length();
}
/**
 *     @param lhs a vector.
 *    @param rhs a vector.
 *    @return lhs > rhs
 */
template<class T,std::size_t Sz>
inline bool operator> (const vector<T,Sz> &lhs, const vector<T,Sz> &rhs){
    return lhs.sqr_length() > rhs.sqr_length();
}
/**
 *    @param lhs a vector.
 *    @param rhs a vector.
 *    @return lhs >= rhs
 */
template<class T,std::size_t Sz>
inline bool operator>= (const vector<T,Sz> &lhs, const vector<T,Sz> &rhs){
    return !(lhs < rhs);
}
/**
 *    @param lhs a vector.
 *    @param rhs a vector.
 *    @return lhs <= rhs
 */

template<class T,std::size_t Sz>
inline bool operator<= (const vector<T,Sz> &lhs, const vector<T,Sz> &rhs){
    return !(lhs > rhs);
}

#endif

    //@}

    //-----------------------------------------------
    //output
    /**
     *    @name output
     *    @relates vector
     */
    //@{

    /**
     *    ostream <<
     */
    template <typename T, std::size_t Sz, typename CharT, class Traits>
    std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const vector<T, Sz>& rhs)
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

    //@}

} //End of namespaces.

#endif
