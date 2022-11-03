// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once


#include "LinearAlgebra.h"
#include "MathExtensions.h"
#include "MultiArrayComplex.h"


namespace quantumdata {


template<typename>
constexpr auto MultiArrayRank_v=std::nullopt;


/// Adding common linear algebra functionality to StateVector and DensityOperator.
template<typename Derived>
struct ArrayBase : linalg::VectorSpace<Derived,cppqedutils::MultiArray<dcomp,MultiArrayRank_v<Derived>>>
{
  using cppqedutils::MultiArray<dcomp,MultiArrayRank_v<Derived>>::MultiArray;
  
  /// \name Naive vector-space operations
  //@{
  Derived& operator+=(const ArrayBase& other) {for (auto&& [t,o] : boost::combine(this->dataView,other.dataView) ) t+=o; return static_cast<Derived&>(*this);}
  Derived& operator-=(const ArrayBase& other) {for (auto&& [t,o] : boost::combine(this->dataView,other.dataView) ) t-=o; return static_cast<Derived&>(*this);}
  //@}

  /// \name Naive vector-space operations allowing also for mixed-mode arithmetics
  //@{
  Derived& operator*=(const auto& dc); // {arrayLow_*=dc; return static_cast<Derived&>(*this);} ///< \tparam OTHER the “other” type in mixed mode

  Derived& operator/=(const auto& dc); // {arrayLow_/=dc; return static_cast<Derived&>(*this);}
  //@}

  /// \name One-dimensional view of the underlying data
  linalg::CVector vectorView(); // what about constness here??? const {return blitzplusplus::unaryArray(arrayLow_);}
  //@}

  /// The entrywise array norm
  /**
   * \f[\norm A _{\text{F}}=\sqrt{\sum_i\,\abs{A_i}^2},\f]
   * with \f$i\f$ running through all the multi-indices.
   */
  double frobeniusNorm() const { return std_ext::ranges::fold( this->dataView | std::views::transform(cppqedutils::sqrAbs) , 0., std::plus{} ); }

};


} // quantumdata


