/// \briefFile{Primitive type wrapper to allow deriving from types that act and feel like primitives. Source: http://github.com/jehugaleahsa/primitive, commit 83d0d1f. Both files http://github.com/jehugaleahsa/primitive/blob/master/arithmetic_traits.h and http://github.com/jehugaleahsa/primitive/blob/master/primitive.hpp are contained in this file. Naming conventions are different here to underline that the code is taken from somewhere else.}

#ifndef CPPQEDCORE_UTILS_PRIMITIVE_HPP_INCLUDED
#define CPPQEDCORE_UTILS_PRIMITIVE_HPP_INCLUDED

#include <iosfwd>

#include <type_traits>

namespace cpputils {

// beginning of file "arithmetic_traits.h"

template<typename TFrom, typename TTo> struct is_promotion : std::false_type {};
template<typename TFrom, typename TTo> struct is_conversion : std::false_type {};

template<> struct is_conversion<int, signed char> : std::true_type {};
template<> struct is_conversion<unsigned int, unsigned char> : std::true_type {};
template<> struct is_conversion<int, short> : std::true_type {};
template<> struct is_conversion<unsigned int, unsigned short> : std::true_type {};

template<> struct is_promotion<signed char, short> : std::true_type {};

template<> struct is_promotion<unsigned char, unsigned short> : std::true_type {};

template<> struct is_promotion<signed char, int> : std::true_type {};
template<> struct is_promotion<short, int> : std::true_type {};

template<> struct is_promotion<unsigned char, unsigned int> : std::true_type {};
template<> struct is_promotion<unsigned short, unsigned int> : std::true_type {};

template<> struct is_promotion<signed char, long> : std::true_type {};
template<> struct is_promotion<short, long> : std::true_type {};
template<> struct is_promotion<int, long> : std::true_type {};

template<> struct is_promotion<unsigned char, unsigned long> : std::true_type {};
template<> struct is_promotion<unsigned short, unsigned long> : std::true_type {};
template<> struct is_promotion<unsigned int, unsigned long> : std::true_type {};

template<> struct is_promotion<signed char, long long> : std::true_type {};
template<> struct is_promotion<short, long long> : std::true_type {};
template<> struct is_promotion<int, long long> : std::true_type {};
template<> struct is_promotion<long, long long> : std::true_type {};

template<> struct is_promotion<unsigned char, unsigned long long> : std::true_type {};
template<> struct is_promotion<unsigned short, unsigned long long> : std::true_type {};
template<> struct is_promotion<unsigned int, unsigned long long> : std::true_type {};
template<> struct is_promotion<unsigned long, unsigned long long> : std::true_type {};

template<> struct is_promotion<signed char, double> : std::true_type {};
template<> struct is_promotion<unsigned char, double> : std::true_type {};
template<> struct is_promotion<short, double> : std::true_type {};
template<> struct is_promotion<unsigned short, double> : std::true_type {};
template<> struct is_promotion<int, double> : std::true_type {};
template<> struct is_promotion<unsigned int, double> : std::true_type {};
template<> struct is_promotion<long, double> : std::true_type {};
template<> struct is_promotion<unsigned long, double> : std::true_type {};
template<> struct is_promotion<float, double> : std::true_type {};

template<> struct is_promotion<signed char, long double> : std::true_type {};
template<> struct is_promotion<unsigned char, long double> : std::true_type {};
template<> struct is_promotion<short, long double> : std::true_type {};
template<> struct is_promotion<unsigned short, long double> : std::true_type {};
template<> struct is_promotion<int, long double> : std::true_type {};
template<> struct is_promotion<unsigned int, long double> : std::true_type {};
template<> struct is_promotion<long, long double> : std::true_type {};
template<> struct is_promotion<unsigned long, long double> : std::true_type {};
template<> struct is_promotion<long long, long double> : std::true_type {};
template<> struct is_promotion<unsigned long long, long double> : std::true_type {};
template<> struct is_promotion<float, long double> : std::true_type {};
template<> struct is_promotion<double, long double> : std::true_type {};

template<> struct is_promotion<signed char, dcomp> : std::true_type {};
template<> struct is_promotion<unsigned char, dcomp> : std::true_type {};
template<> struct is_promotion<short, dcomp> : std::true_type {};
template<> struct is_promotion<unsigned short, dcomp> : std::true_type {};
template<> struct is_promotion<int, dcomp> : std::true_type {};
template<> struct is_promotion<unsigned int, dcomp> : std::true_type {};
template<> struct is_promotion<long, dcomp> : std::true_type {};
template<> struct is_promotion<unsigned long, dcomp> : std::true_type {};
template<> struct is_promotion<float, dcomp> : std::true_type {};
template<> struct is_promotion<double, dcomp> : std::true_type {};


// end of file "arithmetic_traits.h"


template< typename T >
struct IsArithmetic : std::is_arithmetic<T> {};

template<>
struct IsArithmetic<dcomp> : std::true_type {}; // this makes that dcomp is also treated as a primitive arithmetic type


template<typename T, typename = std::enable_if_t< IsArithmetic<T>::value >>
class primitive {
    T m_value;

public:
    using value_type = T;

    constexpr primitive() noexcept: m_value() {}

    template<typename U, typename = std::enable_if_t<
         std::is_same<T, U>::value || is_promotion<U, T>::value
    >>
    constexpr primitive(U const& value) noexcept : m_value(value) {}

    template<typename U, typename = std::enable_if_t< is_promotion<U, T>::value >>
    constexpr primitive(primitive<U> const& other) noexcept : m_value(other.get()) {}

    template<typename U, typename = std::enable_if_t< is_conversion<U, T>::value >>
    constexpr static primitive from(U const& other) noexcept {
        return primitive(T(other));
    }

    primitive(primitive const&) = default;
    primitive(primitive &&) = default;

    primitive& operator=(primitive const&) = default;
    primitive& operator=(primitive &&) = default;

    constexpr T const& get() const noexcept { return m_value; }

    template<typename U = T, typename = std::enable_if_t< !std::is_same<U, bool>::value  >>
    constexpr primitive const& operator+() const noexcept {
        return *this;
    }
    template<typename U = T, typename = std::enable_if_t< !std::is_same<U, bool>::value  >>
    constexpr primitive operator-() const noexcept {
        return primitive(-m_value);
    }

    template<typename U = T, typename = std::enable_if_t< std::is_integral<U>::value && !std::is_same<U, bool>::value >>
    constexpr primitive operator~() const noexcept {
        return primitive(~m_value);
    }

    template<typename U = T, typename = std::enable_if_t< std::is_same<U, bool>::value >>
    constexpr bool operator!() const noexcept {
        return !m_value;
    }

    primitive& operator++() noexcept {
        ++m_value;
        return *this;
    }
    primitive operator++(int) noexcept {
        return primitive(m_value++);
    }

    primitive& operator--() noexcept {
        --m_value;
        return *this;
    }
    primitive operator--(int) noexcept {
        return primitive(m_value--);
    }

    template<typename U>
    primitive& operator+=(U const& other) noexcept {
        m_value += other;
        return *this;
    }
    template<typename U>
    primitive& operator+=(primitive<U> const& other) noexcept {
        m_value += other.get();
        return *this;
    }

    template<typename U>
    primitive& operator-=(U const& other) noexcept {
        m_value -= other;
        return *this;
    }
    template<typename U>
    primitive& operator-=(primitive<U> const& other) noexcept {
        m_value -= other.get();
        return *this;
    }

    template<typename U>
    primitive& operator*=(U const& other) noexcept {
        m_value *= other;
        return *this;
    }
    template<typename U>
    primitive& operator*=(primitive<U> const& other) noexcept {
        m_value *= other.get();
        return *this;
    }

    template<typename U>
    primitive& operator/=(U const& other) noexcept {
        m_value /= other;
        return *this;
    }
    template<typename U>
    primitive& operator/=(primitive<U> const& other) noexcept {
        m_value /= other.get();
        return *this;
    }

    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator%=(U const& other) noexcept {
        m_value %= other;
        return *this;
    }
    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator%=(primitive<U> const& other) noexcept {
        m_value %= other.get();
        return *this;
    }

    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator<<=(U const& other) noexcept {
        m_value <<= other;
        return *this;
    }
    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator<<=(primitive<U> const& other) noexcept {
        m_value <<= other.get();
        return *this;
    }

    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator>>=(U const& other) noexcept {
        m_value >>= other;
        return *this;
    }
    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator>>=(primitive<U> const& other) noexcept {
        m_value >>= other.get();
        return *this;
    }

    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator&=(U const& other) noexcept {
        m_value &= other;
        return *this;
    }
    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator&=(primitive<U> const& other) noexcept {
        m_value &= other.get();
        return *this;
    }

    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator|=(U const& other) noexcept {
        m_value |= other;
        return *this;
    }
    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator|=(primitive<U> const& other) noexcept {
        m_value |= other.get();
        return *this;
    }

    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator^=(U const& other) noexcept {
        m_value ^= other;
        return *this;
    }
    template<typename U, typename = std::enable_if_t< std::is_integral<T>::value && std::is_integral<U>::value >>
    primitive& operator^=(primitive<U> const& other) noexcept {
        m_value ^= other.get();
        return *this;
    }

    template<typename U>
    constexpr explicit operator primitive<U>() const noexcept {
        return primitive<U>(static_cast<U>(m_value));
    }

    friend std::istream& operator>>(std::istream& lhs, primitive<T> & rhs) {
        return lhs >> rhs.m_value;
    }
};

template<typename T>
constexpr primitive<T> operator+(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() + rhs);
}
template<typename T>
constexpr primitive<T> operator+(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs + rhs.get());
}
template<typename T1, typename T2>
constexpr auto operator+(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() + rhs.get())>(lhs.get() + rhs.get());
}

template<typename T>
constexpr primitive<T> operator-(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() - rhs);
}
template<typename T>
constexpr primitive<T> operator-(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs - rhs.get());
}
template<typename T1, typename T2>
constexpr auto operator-(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() - rhs.get())>(lhs.get() - rhs.get());
}

template<typename T>
constexpr primitive<T> operator*(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() * rhs);
}
template<typename T>
constexpr primitive<T> operator*(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs * rhs.get());
}
template<typename T1, typename T2>
constexpr auto operator*(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() * rhs.get())>(lhs.get() * rhs.get());
}

template<typename T>
constexpr primitive<T> operator/(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() / rhs);
}
template<typename T>
constexpr primitive<T> operator/(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs / rhs.get());
}
template<typename T1, typename T2>
constexpr auto operator/(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() / rhs.get())>(lhs.get() / rhs.get());
}

template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator%(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() % rhs);
}
template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator%(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs % rhs.get());
}
template<typename T1, typename T2, typename = std::enable_if_t< std::is_integral<T1>::value && std::is_integral<T2>::value >>
constexpr auto operator%(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() % rhs.get())>(lhs.get() % rhs.get());
}

template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator&(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() & rhs);
}
template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator&(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs & rhs.get());
}
template<typename T1, typename T2, typename = std::enable_if_t< std::is_integral<T1>::value && std::is_integral<T2>::value >>
constexpr auto operator&(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() & rhs.get())>(lhs.get() & rhs.get());
}

template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator|(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() | rhs);
}
template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator|(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs | rhs.get());
}
template<typename T1, typename T2, typename = std::enable_if_t< std::is_integral<T1>::value && std::is_integral<T2>::value >>
constexpr auto operator|(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() | rhs.get())>(lhs.get() | rhs.get());
}

template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator^(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() ^ rhs);
}
template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator^(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs ^ rhs.get());
}
template<typename T1, typename T2, typename = std::enable_if_t< std::is_integral<T1>::value && std::is_integral<T2>::value >>
constexpr auto operator^(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept  {
    return primitive<decltype(lhs.get() ^ rhs.get())>(lhs.get() ^ rhs.get());
}

template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator<<(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() << rhs);
}
template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator<<(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs << rhs.get());
}
template<typename T1, typename T2, typename = std::enable_if_t< std::is_integral<T1>::value && std::is_integral<T2>::value >>
constexpr auto operator<<(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() << rhs.get())>(lhs.get() << rhs.get());
}

template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator>>(primitive<T> const& lhs, T const& rhs) noexcept {
    return primitive<T>(lhs.get() >> rhs);
}
template<typename T, typename = std::enable_if_t< std::is_integral<T>::value >>
constexpr primitive<T> operator>>(T const& lhs, primitive<T> const& rhs) noexcept {
    return primitive<T>(lhs >> rhs.get());
}
template<typename T1, typename T2, typename = std::enable_if_t< std::is_integral<T1>::value && std::is_integral<T2>::value >>
constexpr auto operator>>(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return primitive<decltype(lhs.get() >> rhs.get())>(lhs.get() >> rhs.get());
}

constexpr bool operator&&(primitive<bool> const& lhs, bool const& rhs) noexcept {
    return lhs.get() && rhs;
}
constexpr bool operator&&(bool const& lhs, primitive<bool> const& rhs) noexcept {
    return lhs && rhs.get();
}
constexpr bool operator&&(primitive<bool> const& lhs, primitive<bool> const& rhs) noexcept {
    return lhs.get() && rhs.get();
}

constexpr bool operator||(primitive<bool> const& lhs, bool const& rhs) noexcept {
    return lhs.get() || rhs;
}
constexpr bool operator||(bool const& lhs, primitive<bool> const& rhs) noexcept {
    return lhs || rhs.get();
}
constexpr bool operator||(primitive<bool> const& lhs, primitive<bool> const& rhs) noexcept {
    return lhs.get() || rhs.get();
}

template<typename T>
constexpr bool operator==(primitive<T> const& lhs, T const& rhs) noexcept {
    return lhs.get() == rhs;
}
template<typename T>
constexpr bool operator==(T const& lhs, primitive<T> const& rhs) noexcept {
    return lhs == rhs.get();
}
template<typename T1, typename T2>
constexpr bool operator==(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return lhs.get() == rhs.get();
}

template<typename T>
constexpr bool operator!=(primitive<T> const& lhs, T const& rhs) noexcept {
    return lhs.get() != rhs;
}
template<typename T>
constexpr bool operator!=(T const& lhs, primitive<T> const& rhs) noexcept {
    return lhs != rhs.get();
}
template<typename T1, typename T2>
constexpr bool operator!=(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return lhs.get() != rhs.get();
}

template<typename T>
constexpr bool operator<(primitive<T> const& lhs, T const& rhs) noexcept {
    return lhs.get() < rhs;
}
template<typename T>
constexpr bool operator<(T const& lhs, primitive<T> const& rhs) noexcept {
    return lhs < rhs.get();
}
template<typename T1, typename T2>
constexpr bool operator<(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return lhs.get() < rhs.get();
}

template<typename T>
constexpr bool operator<=(primitive<T> const& lhs, T const& rhs) noexcept {
    return lhs.get() <= rhs;
}
template<typename T>
constexpr bool operator<=(T const& lhs, primitive<T> const& rhs) noexcept {
    return lhs <= rhs.get();
}
template<typename T1, typename T2>
constexpr bool operator<=(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return lhs.get() <= rhs.get();
}

template<typename T>
constexpr bool operator>(primitive<T> const& lhs, T const& rhs) noexcept {
    return lhs.get() > rhs;
}
template<typename T>
constexpr bool operator>(T const& lhs, primitive<T> const& rhs) noexcept {
    return lhs > rhs.get();
}
template<typename T1, typename T2>
constexpr bool operator>(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return lhs.get() > rhs.get();
}

template<typename T>
constexpr bool operator>=(primitive<T> const& lhs, T const& rhs) noexcept {
    return lhs.get() >= rhs;
}
template<typename T>
constexpr bool operator>=(T const& lhs, primitive<T> const& rhs) noexcept {
    return lhs >= rhs.get();
}
template<typename T1, typename T2>
constexpr bool operator>=(primitive<T1> const& lhs, primitive<T2> const& rhs) noexcept {
    return lhs.get() >= rhs.get();
}

template<typename T>
std::ostream& operator<<(std::ostream& lhs, primitive<T> const& rhs) {
    return lhs << rhs.get();
}

} // cpputils

#endif // CPPQEDCORE_UTILS_PRIMITIVE_HPP_INCLUDED
