// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "LinearAlgebra.h"
#include "MultiArrayComplex.h"


namespace quantumdata {


template<typename>
constexpr auto multiArrayRank_v=std::nullopt;


/// Adding common linear algebra functionality to MultiArray for use in StateVector and DensityOperator.
template<typename Derived>
struct ArrayBase : linalg::VectorSpace<Derived,cppqedutils::MultiArray<dcomp,multiArrayRank_v<Derived>>>
{
  ArrayBase(auto&& ... a) : linalg::VectorSpace<Derived,cppqedutils::MultiArray<dcomp,multiArrayRank_v<Derived>>>(std::forward<decltype(a)>(a)...) {}

  /// \name Naive vector-space operations
  //@{
  Derived& operator+=(const ArrayBase& other) {for (auto&& [t,o] : boost::combine(this->dataView,other.dataView) ) t+=o; return static_cast<Derived&>(*this);}
  Derived& operator-=(const ArrayBase& other) {for (auto&& [t,o] : boost::combine(this->dataView,other.dataView) ) t-=o; return static_cast<Derived&>(*this);}
  //@}

  /// \name Naive vector-space operations allowing also for mixed-mode arithmetics
  //@{
  Derived& operator*=(const auto& dc) {for (auto&& t : this->mutableView().dataView) t*=dc; return static_cast<Derived&>(*this);}
  Derived& operator/=(const auto& dc) {for (auto&& t : this->mutableView().dataView) t/=dc; return static_cast<Derived&>(*this);}
  //@}

  /// The entrywise array norm
  /**
   * \f[\norm A _{\text{F}}=\sqrt{\sum_i\,\abs{A_i}^2},\f]
   * with \f$i\f$ running through all the multi-indices.
   */
  friend double frobeniusNorm(const ArrayBase& a) { return sqrt( std_ext::ranges::fold( a.dataView | std::views::transform(sqrAbs) , 0., std::plus{} ) ); }

};


} // quantumdata



/*
namespace boost { namespace numeric { namespace odeint {
  
template<typename Derived>
struct is_resizeable<::quantumdata::ArrayBase<Derived>> : boost::true_type {};

template<typename Derived>
struct same_size_impl<::quantumdata::ArrayBase<Derived>, ::quantumdata::ArrayBase<Derived>>
{ // define how to check size
  static bool same_size(const ::quantumdata::ArrayBase<Derived> &v1, const ::quantumdata::ArrayBase<Derived> &v2) {return v1.extents == v2.extents;}
};

/// TODO: a reserve could be defined for the vector to be resized
template<typename Derived>
struct resize_impl<::quantumdata::ArrayBase<Derived>, ::quantumdata::ArrayBase<Derived>>
{ // define how to resize
  static void resize(::quantumdata::ArrayBase<Derived> &v1, const ::quantumdata::ArrayBase<Derived> &v2) {v1.resize( v2.extents );}
};

template<typename Derived>
struct vector_space_norm_inf<::quantumdata::ArrayBase<Derived>>
{
  typedef double result_type;
  double operator()(const ::quantumdata::ArrayBase<Derived>& v ) const
  {
    return max( abs(v) );
  }
};

template<typename Derived>
struct norm_result_type<::quantumdata::ArrayBase<Derived>> : mpl::identity<double> {};

} } } // boost::numeric::odeint
*/
