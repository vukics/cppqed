// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once


#include "LinearAlgebra.h"
#include "MultiArrayComplex.h"


namespace quantumdata {


template<typename>
constexpr auto multiArrayRank_v=std::nullopt;


/// Adding common linear algebra functionality to StateVector and DensityOperator.
template<typename Derived, template <typename, size_t> typename MA_BASE=cppqedutils::MultiArrayConstView>
struct ArrayBase : linalg::VectorSpace<Derived,cppqedutils::MultiArray<dcomp,multiArrayRank_v<Derived>,MA_BASE>>
{
  using cppqedutils::MultiArray<dcomp,multiArrayRank_v<Derived>,MA_BASE>::MultiArray;
  
  /// \name Naive vector-space operations
  //@{
  Derived& operator+=(const ArrayBase& other) {for (auto&& [t,o] : boost::combine(this->dataView,other.dataView) ) t+=o; return static_cast<Derived&>(*this);}
  Derived& operator-=(const ArrayBase& other) {for (auto&& [t,o] : boost::combine(this->dataView,other.dataView) ) t-=o; return static_cast<Derived&>(*this);}
  //@}

  /// \name Naive vector-space operations allowing also for mixed-mode arithmetics
  //@{
  Derived& operator*=(const auto& dc) {for (auto&& t : this->dataView) t*=dc; return static_cast<Derived&>(*this);}
  Derived& operator/=(const auto& dc) {for (auto&& t : this->dataView) t/=dc; return static_cast<Derived&>(*this);}
  //@}

  /// The entrywise array norm
  /**
   * \f[\norm A _{\text{F}}=\sqrt{\sum_i\,\abs{A_i}^2},\f]
   * with \f$i\f$ running through all the multi-indices.
   */
  double frobeniusNorm() const { return std_ext::ranges::fold( this->dataView | std::views::transform(sqrAbs) , 0., std::plus{} ); }

};


} // quantumdata


