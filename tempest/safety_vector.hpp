#pragma once
#ifndef TM_SAFETY_VECTOR_HPP
#define TM_SAFETY_VECTOR_HPP

/** 
 * @file safety_vector.hpp
 * @brief Defines tempest::safety_vector<T,Sz>.
 *
 * $RCSfile: safety_vector.hpp,v $
 * $Date: 2005/05/09 12:23:20 $
 * $Author: Toru Matsuoka $
 */

#include <cmath>
#include <stdexcept> //
#include <new>       //placement new
#include <iterator>

//#include<ostream>
#include <sstream>

#include "vector_base.hpp"

namespace tempest
{

    namespace safety_vector_util
    {  
        template <class IT>
        void destroy(IT first, IT last)
        {
            typedef typename std::iterator_traits<IT>::value_type T;
            while (first != last)
            {
                first->~T();
                ++first;
            }
        }

        template <class T, std::size_t Sz>
        class array_releaser
        {
        public:
            typedef T* pointer;

        private: //members
            pointer m_first;
            pointer m_constructed;

        public:
            array_releaser() : m_first(static_cast<T*>(operator new(sizeof(T) * Sz))), m_constructed(m_first) {}

            ~array_releaser()
            {
                if (m_first == 0)
                    return;
                destroy(m_first, m_constructed);
                operator delete(m_first);
            }

        private:
            template <class IT>
            void construct_impl(IT a, IT b, pointer first)
            {
                while (a != b)
                {
                    new (static_cast<void*>(m_constructed)) T(*a); //allocator.construct(first,*a);
                    ++a;
                    m_constructed++;
                }
            }

            template <class IT, class TAG> //general
            void construct_tag_dispatch(IT a, IT b, TAG)
            {
                pointer last = m_first + Sz;
                while (m_constructed != last)
                {
                    if (a != b)
                    {
                        new (static_cast<void*>(m_constructed)) T(*a); //allocator.construct(first,*a);
                        ++a;
                    }
                    else
                    {
                        new (static_cast<void*>(m_constructed)) T();
                    }
                    m_constructed++;
                }
            }

            template <class IT> //random access
            void construct_tag_dispatch(IT a, IT b, std::random_access_iterator_tag)
            {
                std::ptrdiff_t dist = b - a;
                if (dist < std::ptrdiff_t(Sz))
                {
                    construct_impl(a, b, m_first); //
                    construct();
                }
                else
                {
                    construct_impl(a, a + Sz, m_first);
                }
            }

        public:
            void construct()
            {
                pointer last = m_first + Sz;
                while (m_constructed != last)
                {
                    new (m_constructed) T();
                    m_constructed++;
                }
            }
            template <class IT>
            void construct(IT a, IT b)
            {
                typedef typename std::iterator_traits<IT>::iterator_category tag;
                construct_tag_dispatch(a, b, tag());
            }

            pointer release()
            {
                pointer p = m_first;
                m_first = m_constructed = 0;
                return p;
            }
        };

        template <class T, std::size_t Sz>
        class calc_vector
        {
        public:
            typedef calc_vector this_type;
            typedef safety_vector_util::array_releaser<T, Sz> array_releaser;

        public:
            calc_vector(const safety_vector<T, Sz>& rhs)
            {
                array_releaser tmp;
                tmp.construct(rhs.begin(), rhs.end());
                element = tmp.release();
            }

            ~calc_vector()
            {
                if (element == 0)
                    return;
                safety_vector_util::destroy(element, element + Sz);
                operator delete(element);
            }

            this_type& neg()
            {
                for (std::size_t i = 0; i < Sz; i++)
                {
                    element[i] *= T(-1);
                }
                return *this;
            }

#define DECLARE_OP_EQUAL_SP(OP, METHOD)                \
    this_type& METHOD(const safety_vector<T, Sz>& rhs) \
    {                                                  \
        for (std::size_t i = 0; i < Sz; i++)           \
        {                                              \
            element[i] OP rhs[i];                      \
        }                                              \
        return *this;                                  \
    }                                                  \
    this_type& METHOD(const this_type& rhs)            \
    {                                                  \
        for (std::size_t i = 0; i < Sz; i++)           \
        {                                              \
            element[i] OP rhs[i];                      \
        }                                              \
        return *this;                                  \
    }

            DECLARE_OP_EQUAL_SP(+=, add)
            DECLARE_OP_EQUAL_SP(-=, sub)
            DECLARE_OP_EQUAL_SP(*=, mul)
            DECLARE_OP_EQUAL_SP(/=, div)

#undef DECLARE_OP_EQUAL_SP

#define DECLARE_OP_EQUAL_SCALER_SP(OP, METHOD) \
    this_type& METHOD(const T& rhs)            \
    {                                          \
        for (std::size_t i = 0; i < Sz; i++)   \
        {                                      \
            element[i] OP rhs;                 \
        }                                      \
        return *this;                          \
    }

            DECLARE_OP_EQUAL_SCALER_SP(*=, mul)
            DECLARE_OP_EQUAL_SCALER_SP(/=, div)

#undef DECLARE_OP_EQUAL_SCALER_SP

            T* release()
            {
                T* p = element;
                element = 0;
                return p;
            }

            void swap(safety_vector<T, Sz>& rhs)
            {
                std::swap(element, rhs.element);
            }

        private:
            T* element;
        };
    }

    /** 
 *  @brief Safety_vector is more safe than vector. 
 *  
 *  @code
 *  tempest::safety_vector<double,3> v;
 *  @endcode
 *
 *  @attention This class is slow.
 */

    template <class T, std::size_t Sz>
    class safety_vector : public vector_base<safety_vector<T, Sz>, T, Sz>
    {
        friend class safety_vector_util::calc_vector<T, Sz>;

    public:
        //-----------------------------------------------
        //type defines
        typedef T value_type;
        typedef safety_vector<T, Sz> this_type;

        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        typedef pointer iterator;             //std::iterator<std::random_access_iterator_tag,T,difference_type>
        typedef const_pointer const_iterator; //std::iterator<std::random_access_iterator_tag,const T,difference_type>
        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

        typedef safety_vector_util::array_releaser<T, Sz> array_releaser;
        typedef safety_vector_util::calc_vector<T, Sz> calc_vector;

        //-----------------------------------------------
    public:                                 //static members...
        static const size_type c_size = Sz; //container size
    public:
        //
        calc_vector calc() const
        {
            return calc_vector(*this);
        }

        this_type clone() const
        {
            return this_type(*this);
        }

        //functions for iterator
        pointer begin() { return element; }
        pointer end() { return element + c_size; }
        const_pointer begin() const { return element; }
        const_pointer end() const { return element + c_size; }

        reverse_iterator rbegin() { return reverse_iterator(end()); }
        reverse_iterator rend() { return reverse_iterator(begin()); }
        const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
        const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }

        //front & back functions
        //
        reference front() { return *begin(); }
        reference back() { return *end(); }
        const_reference front() const { return *begin(); }
        const_reference back() const { return *end(); }

        //constructor
        //
        safety_vector()
        {
            array_releaser tmp;
            tmp.construct();
            element = tmp.release();
        }

        safety_vector(const this_type& rhs)
        {
            array_releaser tmp;
            tmp.construct(rhs.begin(), rhs.end());
            element = tmp.release();
        }

        safety_vector(calc_vector& rhs) : element(0)
        {
            rhs.swap(*this);
        }

        template <class IT>
        safety_vector(IT first, IT last)
        {
            array_releaser tmp;
            tmp.construct(first, last);
            element = tmp.release();
        }

        explicit safety_vector(value_type const(&array)[Sz])
        {
            array_releaser tmp;
            tmp.construct(&(array[0]), &(array[c_size]));
            element = tmp.release();
        }

        template <class Self>
        explicit safety_vector(const vector_base<Self, T, Sz>& rhs)
        {
            array_releaser tmp;
            tmp.construct(static_cast<const Self&>(rhs).begin(), static_cast<const Self&>(rhs).end());
            element = tmp.release();
        }

        ~safety_vector()
        {
            safety_vector_util::destroy(element, element + c_size);
            operator delete(element);
        }

        void swap(this_type& rhs)
        {
            std::swap(element, rhs.element);
        }

        //-----------------------------------------------
        //inserters
        this_type& operator=(const this_type& rhs)
        {
            this_type temp(rhs);
            swap(temp);
            return *this;
        }

        this_type& operator=(calc_vector& rhs)
        {
            rhs.swap(*this);
            return *this;
        }

        template <class Self>
        this_type& operator=(const vector_base<Self, T, Sz>& rhs)
        {
            this_type temp(*this);
            std::copy(static_cast<const Self&>(rhs).begin(), static_cast<const Self&>(rhs).end(), temp.begin());
            swap(temp);
            return *this;
        }

        template <class IT>
        void assign(IT first, IT last)
        {
            swap(this_type(first, last));
        }

        //-----------------------------------------------
        //capacity
        size_type size() const { return Sz; }
        size_type max_size() const { return Sz; }
        bool empty() const { return false; }

        //-----------------------------------------------
        //operators

        this_type& negate()
        {
            this_type temp(*this);
            for (size_type i = 0; i < Sz; i++)
            {
                temp[i] = -temp[i];
            }
            swap(temp);
            return *this;
        }

/**
 *  @def DECLARE_OP_EQUAL( OP )
 *  Defines "safety_vector<T,Sz,Al> OP= safety_vector<T,Sz,Al>".
 */
#define DECLARE_OP_EQUAL(OP)                     \
    this_type& operator OP(const this_type& rhs) \
    {                                            \
        this_type temp(*this);                   \
        for (size_type i = 0; i < Sz; i++)       \
        {                                        \
            temp[i] OP rhs[i];                   \
        }                                        \
        swap(temp);                              \
        return *this;                            \
    }

        DECLARE_OP_EQUAL(+= )
        DECLARE_OP_EQUAL(-= )
        DECLARE_OP_EQUAL(*= )
        DECLARE_OP_EQUAL(/= )

#undef DECLARE_OP_EQUAL

#define DECLARE_OP_EQUAL_SP(OP, METHOD)     \
    this_type& METHOD(const this_type& rhs) \
    {                                       \
        for (size_type i = 0; i < Sz; i++)  \
        {                                   \
            element[i] OP rhs[i];           \
        }                                   \
        return *this;                       \
    }

        DECLARE_OP_EQUAL_SP(+=, add)
        DECLARE_OP_EQUAL_SP(-=, sub)
        DECLARE_OP_EQUAL_SP(*=, mul)
        DECLARE_OP_EQUAL_SP(/=, div)

#undef DECLARE_OP_EQUAL_SP

/**
 *  @def DECLARE_OP_EQUAL_SCALER( OP )
 *  Defines Defines "safety_vector<T,Sz> \OP= T".
 */

#define DECLARE_OP_EQUAL_SCALER(OP)               \
    this_type& operator OP(const value_type& rhs) \
    {                                             \
        this_type temp(*this);                    \
        for (size_type i = 0; i < Sz; i++)        \
        {                                         \
            temp[i] OP rhs;                       \
        }                                         \
        swap(temp);                               \
        return *this;                             \
    }

        DECLARE_OP_EQUAL_SCALER(*= )
        DECLARE_OP_EQUAL_SCALER(/= )

#undef DECLARE_OP_EQUAL_SCALER

#define DECLARE_OP_EQUAL_SCALER_SP(OP, METHOD) \
    this_type& METHOD(const value_type& rhs)   \
    {                                          \
        for (size_type i = 0; i < Sz; i++)     \
        {                                      \
            element[i] OP rhs;                 \
        }                                      \
        return *this;                          \
    }

        DECLARE_OP_EQUAL_SCALER_SP(*=, mul)
        DECLARE_OP_EQUAL_SCALER_SP(/=, div)

#undef DECLARE_OP_EQUAL_SCALER_SP

        //

        T& operator[](size_type i)
        {
            return element[i];
        }

        const T& operator[](size_type i) const
        {
            return element[i];
        }

        T& at(size_type i)
        {
            if (Sz <= i)
            {
                throw std::out_of_range("tempest::safety_vector");
            }
            return element[i];
        }
        const T& at(size_type i) const
        {
            if (Sz <= i)
            {
                throw std::out_of_range("tempest::safety_vector");
            }
            return element[i];
        }

        //-----------------------------------------------
        //utilities

        value_type length() const
        {
            using namespace std;
            return sqrt(sqr_length());
        }

        value_type sqr_length() const
        {
            T temp = T();
            for (size_type i = 0; i < Sz; i++)
            {
                temp += element[i] * element[i];
            }
            return temp;
        }

        value_type sum() const
        {
            T temp = T();
            for (size_type i = 0; i < Sz; i++)
            {
                temp += element[i];
            }
            return temp;
        }

        this_type& normalize()
        {
            using namespace std;

            T length = sqr_length(); //||V||^2
            if (length == T())
                return *this;

            length = T(1) / sqrt(length); // 1 / ||V||
            //length = rsqrt(length);
            swap((this_type(*this) *= length));

            return *this;
        }

        //-----------------------------------------------
        //debug code

    public:
        const char* debug() const { return "tempest::safety_vector<T,Sz>"; }
    private:
        pointer element;
    };

    //-----------------------------------------------
    // utility functions

    template <class T, std::size_t Sz>
    inline T length(const safety_vector<T, Sz>& rhs)
    {
        return rhs.length();
    }

    template <class T, std::size_t Sz>
    inline T sqr_length(const safety_vector<T, Sz>& rhs)
    {
        return rhs.sqr_length();
    }

    template <class T, std::size_t Sz>
    inline T sum(const safety_vector<T, Sz>& rhs)
    {
        return rhs.sum();
    }

    template <class T, std::size_t Sz>
    inline safety_vector<T, Sz> normalize(const safety_vector<T, Sz>& rhs)
    {
        T length = rhs.sqr_length(); //||V||^2
        //if (length == T()) return rhs;

        length = static_cast<T>(1) / std::sqrt(length); // 1 / ||V||
        //length = rsqrt(length);
        return safety_vector<T, Sz>(rhs).mul(length);
    }

    template <class T, std::size_t Sz>
    inline T dot(const safety_vector<T, Sz>& lhs, const safety_vector<T, Sz>& rhs)
    {
        T temp = T();
        for (std::size_t i = 0; i < Sz; i++)
        {
            temp += (lhs[i] * rhs[i]);
        }
        return temp;
    }

    //-----------------------------------------------
    //unary operators
    template <class T, std::size_t Sz>
    inline safety_vector<T, Sz> operator+(const safety_vector<T, Sz>& rhs)
    {
        return rhs;
    }
    template <class T, std::size_t Sz>
    inline safety_vector<T, Sz> operator-(const safety_vector<T, Sz>& rhs)
    {
        return safety_vector<T, Sz>(rhs).mul(T(-1));
    }

//-----------------------------------------------
//binary operators

#define DECLARE_OP_BIN(OP, METHOD)                                                \
    \
template<class T, std::size_t Sz> \
inline safety_vector<T, Sz>                                                       \
    operator OP(const safety_vector<T, Sz>& lhs, const safety_vector<T, Sz>& rhs) \
    {                                                                             \
        return safety_vector<T, Sz>(lhs).METHOD(rhs);                             \
    \
}

    DECLARE_OP_BIN(+, add)
    DECLARE_OP_BIN(-, sub)
    DECLARE_OP_BIN(*, mul)
    DECLARE_OP_BIN(/, div)

#undef DECLARE_OP_BIN

    /*
template<class T,std::size_t Sz>    
inline calc_vector operator+(const     safety_vector<T,Sz> &lhs, const     safety_vector<T,Sz> &rhs){//sf sf
    return lhs.swapper().add(rhs);
}
    
template<class T,std::size_t Sz>
inline calc_vector operator + (const     safety_vector<T,Sz> &lhs, const     calc_vector &rhs){//sf sw
    return calc_vector(lhs) . add (rhs);
}
    
template<class T,std::size_t Sz>
inline calc_vector& operator + (calc_vector &lhs,    const     safety_vector<T,Sz> &rhs){//sw sf
    return lhs . add (rhs);
}
    
template<class T,std::size_t Sz>
inline calc_vector& operator + ( calc_vector &lhs,    const     calc_vector &rhs){//sw sw
    return lhs . add (rhs);
}
*/

    //-----------------------------------------------
    //specific scalar

    template <class T, std::size_t Sz>
    inline safety_vector<T, Sz> operator*(const T lhs, const safety_vector<T, Sz>& rhs)
    {
        return safety_vector<T, Sz>(rhs).mul(lhs);
    }

    template <class T, std::size_t Sz>
    inline safety_vector<T, Sz> operator*(const safety_vector<T, Sz>& lhs, const T rhs)
    {
        return safety_vector<T, Sz>(lhs).mul(rhs);
    }

    template <class T, std::size_t Sz>
    inline safety_vector<T, Sz> operator/(const safety_vector<T, Sz>& lhs, const T rhs)
    {
        return safety_vector<T, Sz>(lhs).div(rhs);
    }

    //--------------------------------------------------
    //compare

    template <class T, std::size_t Sz>
    inline bool operator==(const safety_vector<T, Sz>& lhs, const safety_vector<T, Sz>& rhs)
    {
        for (std::size_t i = 0; i < Sz; ++i)
        {
            if (lhs[i] != rhs[i])
                return false;
        }
        return true;
    }
    template <class T, std::size_t Sz>
    inline bool operator!=(const safety_vector<T, Sz>& lhs, const safety_vector<T, Sz>& rhs)
    {
        return !(lhs == rhs);
    }

#if 0
template<class T,std::size_t Sz,class Al>
inline bool operator< (const safety_vector<T,Sz,Al> &lhs,const safety_vector<T,Sz,Al> &rhs){
    return lhs.sqr_length() < rhs.sqr_length();
}
template<class T,std::size_t Sz,class Al>
inline bool operator> (const safety_vector<T,Sz,Al> &lhs,const safety_vector<T,Sz,Al> &rhs){
    return lhs.sqr_length() > rhs.sqr_length();
}
template<class T,std::size_t Sz,class Al>
inline bool operator<= (const safety_vector<T,Sz,Al> &lhs,const safety_vector<T,Sz,Al> &rhs){
    return !(lhs > rhs);
}
template<class T,std::size_t Sz,class Al>
inline bool operator>= (const safety_vector<T,Sz,Al> &lhs,const safety_vector<T,Sz,Al> &rhs){
    return !(lhs < rhs);
}

#endif

    //-----------------------------------------------
    //output

    /** 
 * ostream << 
 */
    template <typename T, std::size_t Sz, typename _CharT, class _Traits>
    std::basic_ostream<_CharT, _Traits>& operator<<(std::basic_ostream<_CharT, _Traits>& os, const safety_vector<T, Sz>& rhs)
    {

        std::basic_ostringstream<_CharT, _Traits> s;
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
}
#endif
