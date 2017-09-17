#pragma once
#ifndef TM_ET_VECTOR_HPP
#define TM_ET_VECTOR_HPP

/**
 *  @file et_vector.hpp
 *  @author Toru Matsuoka
 *  @brief Defines structs and functions related for ET(Expression Template).
 *
 *  ET(Expression Template) demands from compiler high optimizations.
 *  But, in fact, a fine compiler would generate same code whether you used ET
 * or not.
 *  I think ... ET is wonderful hack, but it is only for code guru.
 *  @see http://osl.iu.edu/~tveldhui/papers/Expression-Templates/exprtmpl.html
 *
 * $Id: et_vector.hpp,v 1.5 2005/02/26 17:16:53      Exp $
 */
//

#include "vector_base.hpp"
#include "operation_result.hpp"
#include <cstddef>
#include <ostream>
#include <sstream>
#include <limits>

#define MAX_SCALAR_VECTOR_SIZE 255
//#define MAX_SCALAR_VECTOR_SIZE ((std::size_t)-2)

#define MEMBER(TYPE,VAR) const TYPE VAR
//#define INIT(TYPE,VAR)   const TYPE VAR
//#define CALL(TYPE,VAR)   const TYPE VAR

//#define MEMBER(TYPE, VAR) const TYPE& VAR
#define INIT(TYPE, VAR) const TYPE& VAR
#define CALL(TYPE, VAR) const TYPE& VAR

namespace tempest
{
    namespace ns_et
    {
        struct posite
        {
            template <class T>
            static inline const T apply(const T& a)
            {
                return a;
            }
        };
        struct negate
        {
            template <class T>
            static inline const T apply(const T& a)
            {
                return -a;
            }
        };

/**
 *  @def DECLARE_OP_STRUCT
 *  Defines Expression Template evaluater struct.
 */

#define DECLARE_OP_STRUCT(STRUCT_NAME, OP)                                                \
    struct STRUCT_NAME                                                                    \
    {                                                                                     \
        template <class P, class Q>                                                       \
        static inline typename operation_result<P, Q>::type apply(const P& a, const Q& b) \
        {                                                                                 \
            return a OP b;                                                                \
        }                                                                                 \
    };

        DECLARE_OP_STRUCT(plus, +)
        DECLARE_OP_STRUCT(minus, -)
        DECLARE_OP_STRUCT(multiplies, *)
        DECLARE_OP_STRUCT(divides, / )
        DECLARE_OP_STRUCT(modulus, % )

#undef DECLARE_OP_STRUCT

        template <class Self>
        struct expression_base
        {
            typedef Self this_type;
        };

        template <class Self, class T, std::size_t Sz>
        struct vector_expression_base : expression_base<Self>
        {
            typedef T value_type;
            value_type operator[](std::size_t i) const { return (static_cast<const Self&>(*this))[i]; }
            size_t size()const{return (static_cast<const Self&>(*this)).size();}
        };

        /**
         * @class scalar_vector
         * @brief scalar_vector will be used for Expression Template Node.
         */
        template <class T, std::size_t Sz = MAX_SCALAR_VECTOR_SIZE>
        struct scalar_vector : public vector_expression_base< scalar_vector<T, Sz>, T, Sz >
        {
            // type define
            typedef T value_type;
            // size
            static const std::size_t c_size = Sz;

            explicit scalar_vector(INIT(T, rhs)) : m(rhs) { }

            value_type operator[](std::size_t i) const { return m; }
            size_t size()const{return c_size;}
            // member...
            MEMBER(T, m);
        };

        /**
         * @class unop_vector
         * @brief unop_vector will be used for Expression Template Node.
         */
        template <class OP, class V>
        struct unop_vector : public vector_expression_base<unop_vector<OP, V>, typename V::value_type, V::c_size >
        {
            // type define
            typedef typename V::value_type value_type;
            // size
            static const std::size_t c_size = V::c_size;

            explicit unop_vector(INIT(V, rhs)) : m(rhs){};

            value_type operator[](std::size_t i) const { return OP::apply(m[i]); }
            size_t size()const{return c_size;}
            // member...
            MEMBER(V, m);
        };

        /**
         * @class binop_vector
         * @brief binop_vector will be used for Expression Template Node.
         */
        template <class OP, class P, class Q>
        struct binop_vector : 
            public 
            vector_expression_base<
                binop_vector<OP, P, Q>,
                typename operation_result<
                    typename P::value_type,
                    typename Q::value_type
                >::type,
                (P::c_size < Q::c_size) ? P::c_size : Q::c_size
            >
        {
            typedef typename P::value_type PV;
            typedef typename Q::value_type QV;
            typedef typename operation_result<PV, QV>::type value_type;

            static const std::size_t p_size = P::c_size;
            static const std::size_t q_size = Q::c_size;
            static const std::size_t c_size = (p_size < q_size) ? p_size : q_size;

            binop_vector(INIT(P, lhs), INIT(Q, rhs)) : left(lhs), right(rhs) {}

            value_type operator[](std::size_t i) const { return OP::apply(left[i], right[i]); }
            size_t size()const{return c_size;}
            // member
            MEMBER(P, left);
            MEMBER(Q, right);
        };
        //-----------------------------------------------

        template <class XX>
        inline const expression_base<XX>& operator+(const expression_base<XX>& rhs)
        {
            return rhs;
        }

        template <class XX>
        inline unop_vector<negate, XX> operator-(const expression_base<XX>& rhs)
        {
            return unop_vector<negate, XX>(static_cast<const XX&>(rhs));
        }

#define DECLARE_ET_VECTOR_OPERATOR(STRUCT_NAME, OP)                                                                      \
    template <class TT, class XX, class T, std::size_t SzT, std::size_t SzX>                                             \
    inline binop_vector<STRUCT_NAME, TT, XX> operator OP(const vector_expression_base<TT, T, SzT>& lhs, const vector_expression_base<XX, T, SzX>& rhs) \
    {                                                                                                                    \
        return binop_vector<STRUCT_NAME, TT, XX>(static_cast<const TT&>(lhs), static_cast<const XX&>(rhs));              \
    }                                                                                           

        DECLARE_ET_VECTOR_OPERATOR(plus, +)
        DECLARE_ET_VECTOR_OPERATOR(minus, -)
        DECLARE_ET_VECTOR_OPERATOR(multiplies, *)
        DECLARE_ET_VECTOR_OPERATOR(divides, / )
        //DECLARE_ET_VECTOR_OPERATOR(modulus, % )

#undef DECLARE_ET_VECTOR_OPERATOR

#define DECLARE_ET_VECTOR_OPERATOR(STRUCT_NAME, OP)\
    template <class Self, class T, size_t Sz>\
    inline binop_vector<STRUCT_NAME, Self, scalar_vector<T, Sz> > operator OP(const vector_expression_base<Self, T, Sz>& lhs, const T rhs) \
    {                                                                                                                                      \
        return binop_vector<STRUCT_NAME, Self, scalar_vector<T, Sz> >(static_cast<const Self&>(lhs), scalar_vector<T, Sz>(rhs));           \
    }\
    template <class Self, class T, size_t Sz>\
    inline binop_vector<STRUCT_NAME, scalar_vector<T, Sz>, Self > operator OP(const T lhs, const vector_expression_base<Self, T, Sz>& rhs) \
    {                                                                                                                                      \
        return binop_vector<STRUCT_NAME, scalar_vector<T, Sz>, Self >(scalar_vector<T, Sz>(lhs), static_cast<const Self&>(rhs) );          \
    }          

        DECLARE_ET_VECTOR_OPERATOR(multiplies, *)
        DECLARE_ET_VECTOR_OPERATOR(divides, / )

#undef DECLARE_ET_VECTOR_OPERATOR               

    } // end of namespace et

    //-------------------------------------------------------------------------------------------------
    // Helper for making ref_vector.
    template <class Self, class T, std::size_t Sz>
    inline ns_et::unop_vector<ns_et::posite, Self> et(const vector_base<Self, T, Sz>& rhs)
    {
        return ns_et::unop_vector<ns_et::posite, Self>(static_cast<const Self&>(rhs));
    }

#define DECLARE_SCALAR_ET(TYPE)                    \
    inline ns_et::scalar_vector<TYPE> et(TYPE rhs) \
    {                                              \
        return ns_et::scalar_vector<TYPE>(rhs);    \
    }

    // typedef float TMFLOAT;

    DECLARE_SCALAR_ET(signed char)
    DECLARE_SCALAR_ET(signed short)
    DECLARE_SCALAR_ET(signed int)
    DECLARE_SCALAR_ET(signed long int)
    DECLARE_SCALAR_ET(signed long long int)

    DECLARE_SCALAR_ET(unsigned char)
    DECLARE_SCALAR_ET(unsigned short)
    DECLARE_SCALAR_ET(unsigned int)
    DECLARE_SCALAR_ET(unsigned long int)
    DECLARE_SCALAR_ET(unsigned long long int)

    DECLARE_SCALAR_ET(float)
    DECLARE_SCALAR_ET(double)
    DECLARE_SCALAR_ET(long double)

#undef DECLARE_SCALAR_ET

    template <class _CharT, class _Traits, class SelfE, class T, std::size_t Sz>
    std::basic_ostream<_CharT, _Traits>& operator<<(std::basic_ostream<_CharT, _Traits>& os, const ns_et::vector_expression_base<SelfE, T, Sz>& rhs)
    {
        const size_t sz = Sz;
        std::basic_ostringstream<_CharT, _Traits> s;
        s.flags(os.flags());
        s.imbue(os.getloc());
        s.precision(os.precision());
        s << "(";
        for (std::size_t i = 0; i < sz - 1; ++i)
        {
            s << rhs[i] << ",";
        }
        s << rhs[sz - 1];
        s << ")";
        return os << s.str();
    }

    /*
     * @brief Define  "vector << [vector_expression]".
     * Why I don't define "vector = [vector_expression]" ?
     * This ans is that "operator=" must be non static member function !
     */
    template <class Self, class SelfE, class T, std::size_t Sz>
    inline vector_base<Self, T, Sz>& operator<<(vector_base<Self, T, Sz>& out, const ns_et::vector_expression_base<SelfE, T, Sz>& exp)
    {
        Self& out_ = static_cast<Self&>(out);
        const SelfE& exp_ = static_cast<const SelfE&>(exp);
        const size_t sz = Sz;
        for (size_t i = 0; i < sz; i++)
        {
            out_[i] = exp_[i];
        }
        return out;
    }

} // end of namespace tempest

#undef MEMBER
#undef INIT
#undef CALL

#endif
