// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include <array>
#include <stdexcept>


/// Thrown in the case of dimensionality mismatch of constructs over the same Hilbert space
struct DimensionalityMismatchException : std::invalid_argument {using std::invalid_argument::invalid_argument;};


/// Stores and manipulates dimensions of constructs over composite Hilbert spaces of arbitrary arity
/**
 * \tparamRANK
 */
template<size_t RANK>
class DimensionsBookkeeper
{
public:
  static constexpr size_t N_RANK=RANK; ///< Arity of the Hilbert space

  using Dimensions = std::array<size_t,RANK>; ///< The dimensions as a static vector of size N_RANK

  explicit DimensionsBookkeeper(const Dimensions& dimensions) : dimensions_{dimensions}, totalDimension_{product(dimensions)} {}

  const Dimensions& getDimensions() const {return dimensions_;} ///< Get the Dimensions vector
  
  size_t getTotalDimension() const {return totalDimension_;} ///< Get the total dimension of a system of arbitrary arity

  size_t getDimension() requires ( RANK==1 ) const {return totalDimension_;} ///< Get the (single) dimension for a unary system

  size_t getDimension(size_t i) const {return dimensions_[i];}

  void setDimensions(const Dimensions& dimensions) {dimensions_=dimensions; totalDimension_=product(dimensions);}

private:
  Dimensions dimensions_;
  size_t totalDimension_ ;

};


/// dimensionality comparison for types derived from DimensionsBookkeeper \related DimensionsBookkeeper
template<size_t RANK>
inline bool
operator==(const DimensionsBookkeeper<RANK>& d1, const DimensionsBookkeeper<RANK>& d2)
{
  return d1.getDimensions()==d2.getDimensions();
}

/// dimensionality comparison for types derived from DimensionsBookkeeper \related DimensionsBookkeeper
template<size_t RANK>
inline bool
operator!=(const DimensionsBookkeeper<RANK>& d1, const DimensionsBookkeeper<RANK>& d2)
{
  return !(d1==d2);
}


