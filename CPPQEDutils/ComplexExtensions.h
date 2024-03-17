// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Additional helpers for dcomp}
#pragma once

#include <boost/json.hpp>

#include <boost/serialization/complex.hpp>

#include <complex>


/// Double-precision complex number
/** Even though it is a type, we name it this way because we would like it to closely resemble built-in types */
typedef std::complex<double> dcomp;

using std::complex_literals::operator""i;

inline bool isNonZero(dcomp c) {return bool(real(c)) || bool(imag(c));}

inline bool hasRealPart(dcomp c) {return bool(real(c));}
inline bool hasImagPart(dcomp c) {return bool(imag(c));}

inline bool  absCompare(dcomp c1, dcomp c2) {return  abs(c1)< abs(c2);}
inline bool realCompare(dcomp c1, dcomp c2) {return real(c1)<real(c2);}

/// Generic sqr from this Q&A: https://stackoverflow.com/a/64849248/1171157
auto sqr(auto&& x) // return type is non-reference or trailing
noexcept(noexcept(x*x)) // propagate noexcept
-> decltype(x*x) // enable SFINAE
{ return x * x; }


double relativeDeviation(const auto& a, const auto& b) {return std::abs(a-b)/(std::abs(a)+std::abs(b));}

inline double sqrAbs(dcomp v) {return sqr(std::abs(v));}


namespace boost::json {
    template<class T>
    void tag_invoke(value_from_tag, value& jv, const std::complex<T>& c) {
        jv = {c.real(), c.imag()};
    }

    template<class T>
    std::complex<T> tag_invoke(value_to_tag<std::complex<T>>, const value& jv) {
        return {jv.at(0), jv.at(1)};
    }
} // boost::json


// namespace std {
//
// inline void to_json( ::cppqedutils::json& jv, const complex<double>& dc ) { jv = ::cppqedutils::json{dc.real(),dc.imag()}; }
//
// inline void from_json( const ::cppqedutils::json& jv, complex<double>& dc ) { dc.real(jv[0]); dc.imag(jv[1]); }
//
// }
