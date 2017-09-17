#pragma once
#ifndef TM_META_HPP
#define TM_META_HPP

namespace tempest
{

    /**
     * @struct if_c
     * This struct supports meta programming as "if".
     */
    template <bool b, class P, class Q>
    struct if_c
    {
        typedef P type;
    };
    
    template <class P, class Q>
    struct if_c<false, P, Q>
    {
        typedef Q type;
    };

    //-----------------------------------------------

    //equal boost::is_same<T,U>
    template <class P, class Q>
    struct is_same
    {
        static const bool value = false;
    };
    template <class P>
    struct is_same<P, P>
    {
        static const bool value = true;
    };

    //if(  (T==P)||(T==Q) )return true;
    //else return false;
    template <class T, class P, class Q>
    struct is_contained
    {
        static const bool value = (is_same<T, P>::value | is_same<T, Q>::value);
    };

    //----------------------------------------------

    template <class T>
    struct is_float
    {
        static const bool value = false;
    };

    template <>
    struct is_float<float>
    {
        static const bool value = true;
    };

    template <>
    struct is_float<double>
    {
        static const bool value = true;
    };

    template <>
    struct is_float<long double>
    {
        static const bool value = true;
    };

    template <class T>
    struct is_arith
    {
        static const bool value = std::numeric_limits<T>::is_integer | is_float<T>::value;
    };

    template <class T, bool is_param = false>
    struct param_type_imp
    {
        typedef const T& type;
    };

    template <class T>
    struct param_type_imp<T, true>
    {
        typedef T type;
    };

    template <class T>
    struct param_type
    {
        typedef typename param_type_imp<T, is_arith<T>::value>::type type;
    };
}

#endif
