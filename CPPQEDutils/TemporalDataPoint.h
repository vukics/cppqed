// Copyright András Vukics 2006–2024. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Traits.h"

#include <iostream>

#include <complex>

/// Double-precision complex number
/** Even though it is a type, we name it this way because we would like it to closely resemble built-in types */
typedef std::complex<double> dcomp;


namespace cppqedutils {


// statically_labelled_value
namespace slv_ns {

// fixed_string after https://vector-of-bool.github.io/2021/10/22/string-templates.html
template <size_t Length>
struct fs {
  char chars[Length+1] = {}; // +1 for null terminator

  // Constructor to initialize the _chars array
  constexpr fs(const char (&arr)[Length+1])
  {
    for (size_t i = 0; i < Length; ++i) chars[i] = arr[i];
    chars[Length] = '\0'; // Ensure null termination
  }

};

template <size_t N>
fs(const char (&arr)[N]) -> fs<N-1>;  // Drop the null terminator

template <fs L, typename T>
struct _
{
  _() = default;
  _(auto t) : value(t) {}

  operator T&() {return value;}
  operator const T&() const {return value;}

  T value;

  static constexpr std::string label=L;
};

} // slv_ns


template<slv_ns::fs L, typename T>
auto slv(T t) {return slv_ns::_<L,T>(t);}


/// hana::tuple of double/dcomp or such hana::tuples (recursive definition)
/** Workaround for defining concept recursively from [this Q&A](https://stackoverflow.com/questions/56741456/how-to-define-a-recursive-concept) */
namespace traits {

template <typename T> constexpr bool tdp = false;

template <slv_ns::fs label> constexpr bool tdp<slv_ns::_<label,double>> = true;

template <slv_ns::fs label> constexpr bool tdp<slv_ns::_<label,dcomp>> = true;

template <hana_sequence S> constexpr bool tdp<S> = !!hana::all_of(
  decltype(hana::transform(std::declval<S>(), hana::typeid_)){},
  []<class T>(T) { return tdp<typename T::type>; });

} // traits

template <typename T> concept temporal_data_point = traits::tdp<T>;


static_assert(temporal_data_point<slv_ns::_<"saturation",double>>);
static_assert(temporal_data_point<slv_ns::_<"polarization",dcomp>>);


struct plusTDP
{
  template <slv_ns::fs label>
  slv_ns::_<label,double> operator()(slv_ns::_<label,double> a, slv_ns::_<label,double> b) const {return a.value+b.value;}

  template <slv_ns::fs label>
  slv_ns::_<label,dcomp> operator()(slv_ns::_<label,dcomp> a, slv_ns::_<label,dcomp> b) const {return a.value+b.value;}

  template <temporal_data_point T>
  T operator()(const T& a, const T& b) const {
    using namespace hana::literals;
    T res;
    hana::for_each(hana::range_c<int,0,hana::size(res)>, [&,this] (auto i) {res[i]=(*this)(a[i],b[i]);});
    return res;
  }

};


template <slv_ns::fs label>
std::ostream& streamTDP(slv_ns::_<label,double> tdp, std::ostream& os) {return os<<tdp.value;}

template <slv_ns::fs label>
std::ostream& streamTDP(slv_ns::_<label,dcomp > tdp, std::ostream& os) {return os<<tdp.value;}

inline std::ostream& streamTDP(hana::tuple<>, std::ostream& os) {return os;}

std::ostream& streamTDP(const temporal_data_point auto& tdp, std::ostream& os)
{
  size_t n{0};
  hana::for_each( tdp, [&] (const auto& v) { streamTDP(v,os) << (++n != hana::size(tdp) ? "\t" : "") ; } );
  return os;
}


template <slv_ns::fs label>
void renormTDP(slv_ns::_<label,double>& tdp, double norm) {tdp.value/=norm;}

template <slv_ns::fs label>
void renormTDP(slv_ns::_<label,dcomp >& tdp, double norm) {tdp.value/=norm;}

inline void renormTDP(hana::tuple<>& , double ) {}

void renormTDP(temporal_data_point auto& tdp, double norm) { hana::for_each( tdp, [=] (auto& v) { renormTDP(v,norm) ; } ); }


template <slv_ns::fs label, typename T>
json::value jsonizeTDP_labels(slv_ns::_<label,T>) {return label.chars;}

json::array jsonizeTDP_labels(const temporal_data_point auto& tdp)
{
  json::array res;
  hana::for_each( tdp, [&] (auto v) { res.push_back(jsonizeTDP_labels(v)); } );
  return res;
}

} // cppqedutils
