#pragma once
#ifndef TM_RATIONAL_HPP
#define TM_RATIONAL_HPP

namespace tempest
{

    template <class T>
    class rational
    {
    public:
        rational() {}
        rational(T a, T b = T(1)) : a_(a), b_(b) {}
        rational(const rational& rhs) : a_(rhs.a_), b_(rhs.b_) {}

        rational& operator=(const rational& rhs)
        {
            a_ = rhs.a_;
            b_ = rhs.b_;
            return *this;
        }

        rational& operator+=(const rational& rhs)
        {
            a_ = a_ * rhs.b_ + b_ * rha.a_;
            b_ *= rhs.b_;
            return *this;
        }

        rational& operator-=(const rational& rhs)
        {
            a_ = a_ * rhs.b_ - b_ * rha.a_;
            b_ *= rhs.b_;
            return *this;
        }

        rational& operator*=(const rational& rhs)
        {
            a_ *= rhs.a_;
            b_ *= rhs.b_;
            return *this;
        }

        rational& operator/=(const rational& rhs)
        {
            a_ *= rhs.b_;
            b_ *= rhs.a_;
            return *this;
        }

        rational& operator+=(const T& rhs)
        {
            a_ += rhs * b_;
            return *this;
        }

        rational& operator-=(const T& rhs)
        {
            a_ -= rhs * b_;
            return *this;
        }

        rational& operator*=(const T& rhs)
        {
            a_ *= rhs;
            return *this;
        }

        rational& operator/=(const T& rhs)
        {
            b_ *= rhs;
            return *this;
        }

        operator T() const { return a_ / b_; }

        bool equal(const rational& rhs)
        {
            return (a_ * rhs.b_ == rhs.a_ * b_);
        }

        bool less(const rational& rhs)
        {
            if (b_ * rhs.b_ < 0)
            {
                return a_ * rhs.b_ > rhs.a_ * b_;
            }
            else
            {
                return a_ * rhs.b_ < rhs.a_ * b_;
            }
        }

    private:
        T a_;
        T b_;
    };

#define DECLARE_OPERATOR(OP)                                     \
    \
template<class T> \
inline rational                                                  \
    operator##OP(const rational<T>& lhs, const rational<T>& rhs) \
    {                                                            \
        return rational<T>(lhs) OP## = rhs;                      \
    \
}                                                         \
    \
template<class T> \
inline rational                                                  \
    operator##OP(const rational<T>& lhs, const T& rhs)           \
    {                                                            \
        return rational<T>(lhs) OP## = rhs;                      \
    \
}

    DECLARE_OPERATOR(+)
    DECLARE_OPERATOR(-)
    DECLARE_OPERATOR(*)
    DECLARE_OPERATOR(/ )

#undef DECLARE_OPERATOR

    template <class T>
    inline bool operator==(const rational<T>& lhs, const rational<T>& rhs)
    {
        return lhs.equal(rhs);
    }
    template <class T>
    inline bool operator!=(const rational<T>& lhs, const rational<T>& rhs)
    {
        return !(lhs == rhs);
    }
    template <class T>
    inline bool operator<(const rational<T>& lhs, const rational<T>& rhs)
    {
        return lhs.less(rhs);
    }
    template <class T>
    inline bool operator>=(const rational<T>& lhs, const rational<T>& rhs)
    {
        return !(lhs < rhs);
    }
    template <class T>
    inline bool operator>(const rational<T>& lhs, const rational<T>& rhs)
    {
        return (lhs != rhs) && (lhs >= rhs);
    }
    template <class T>
    inline bool operator<=(const rational<T>& lhs, const rational<T>& rhs)
    {
        return !(lhs > rhs);
    }
}

#endif