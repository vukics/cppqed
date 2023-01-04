// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "MultiArray.h"

#include <stdexcept>


/// Thrown in the case of dimensionality mismatch of constructs over the same Hilbert space
struct DimensionalityMismatchException : std::invalid_argument {using std::invalid_argument::invalid_argument;};


/// Stores and manipulates dimensions of constructs over composite Hilbert spaces of arbitrary arity
template<size_t RANK>
class DimensionsBookkeeper
{
public:
  static constexpr size_t N_RANK=RANK; ///< Arity of the Hilbert space

  using Dimensions = cppqedutils::Extents<RANK>;

  auto operator<=>(const DimensionsBookkeeper&) const = default;
  
  explicit DimensionsBookkeeper(const Dimensions& dimensions) : dimensions_{dimensions} {}

  const Dimensions& getDimensions() const {return dimensions_;} ///< Get the Dimensions vector
  
  size_t getTotalDimension() const {return cppqedutils::multiarray::calculateExtent(dimensions_);} ///< Get the total dimension of a system of arbitrary arity

  size_t getDimension() requires ( RANK==1 ) const {return dimensions_[0];} ///< Get the (single) dimension for a unary system

  size_t getDimension(size_t i) const {return dimensions_[i];}

  void setDimensions(const Dimensions& dimensions) {dimensions_=dimensions;}

private:
  Dimensions dimensions_;

};
