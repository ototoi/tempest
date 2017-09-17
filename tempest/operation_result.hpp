#pragma once
#ifndef TM_OPERATION_RESULT_HPP
#define TM_OPERATION_RESULT_HPP

#include <cstddef>
#include <limits>
#include "meta.hpp"

namespace tempest
{
    template <class T>
    struct type_promotion
    {
        typedef typename if_c<
            (sizeof(T) > sizeof(int)),
            T,
            typename if_c<(std::numeric_limits<T>::digits <= std::numeric_limits<int>::digits),
                          int,
                          unsigned int>::type>::type type;
    };

    template <>
    struct type_promotion<bool>
    {
        typedef int type;
    };

    template <>
    struct type_promotion<wchar_t>
    {
        typedef if_c<(std::numeric_limits<wchar_t>::digits <= std::numeric_limits<int>::digits),
                     int,
                     if_c<(std::numeric_limits<wchar_t>::digits <= std::numeric_limits<unsigned int>::digits),
                          unsigned int,
                          if_c<(std::numeric_limits<wchar_t>::digits <= std::numeric_limits<long int>::digits),
                               long int,
                               unsigned long>::type>::type>::type type;
    };

    template <class P, class Q>
    struct integer_ope_result
    {
        typedef typename if_c<
            is_contained<unsigned long int, P, Q>::value,

            unsigned long int,
            typename if_c<
                ((is_same<long int, P>::value) & (is_same<unsigned int, Q>::value)) | ((is_same<unsigned int, P>::value) & (is_same<long int, Q>::value)),

                typename if_c<
                    (std::numeric_limits<unsigned int>::digits <= std::numeric_limits<long int>::digits),

                    long int,
                    unsigned long int>::type,
                typename if_c<
                    is_contained<long, P, Q>::value,

                    long,
                    typename if_c<
                        is_contained<unsigned, P, Q>::value,

                        unsigned,
                        int>::type>::type>::type>::type type;
    };

    //
    template <class P, class Q>
    struct operation_result
    {
        typedef typename if_c<
            is_contained<long double, P, Q>::value,
            long double,
            typename if_c<
                is_contained<double, P, Q>::value,
                double,
                typename if_c<
                    is_contained<float, P, Q>::value,
                    float,
                    typename integer_ope_result<typename type_promotion<P>::type, typename type_promotion<Q>::type>::type>::type>::type>::type type;
    };

    template <class P, class Q>
    struct return_type2
    {
        typedef typename operation_result<P, Q>::type type;
    };
}

#endif
