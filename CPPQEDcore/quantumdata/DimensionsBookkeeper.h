// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED

#include "BlitzTiny.h"
#include "Exception.h"

#include <type_traits>

/// Thrown in the case of dimensionality mismatch of constructs over the same Hilbert space
class DimensionalityMismatchException : cpputils::Exception {};


/// Stores and manipulates dimensions of constructs over composite Hilbert spaces of arbitrary arity
/**
 * \tparamRANK
 */
template<int RANK>
class DimensionsBookkeeper
{
public:
  static const int N_RANK=RANK; ///< Arity of the Hilbert space

  typedef ExtTiny<RANK> Dimensions; ///< The dimensions as a static vector of size N_RANK

  DimensionsBookkeeper() : dimensions_(), totalDimension_() {}

  explicit DimensionsBookkeeper(const Dimensions& dimensions) : dimensions_(dimensions), totalDimension_(product(dimensions)) {} ///< Straightforward constructor

  const Dimensions& getDimensions    () const {return     dimensions_;} ///< Get the Dimensions vector
  size_t            getTotalDimension() const {return totalDimension_;} ///< Get the total dimension of a system of arbitrary arity

  template<int R=RANK>
  std::enable_if_t<R==1,size_t> getDimension() const {return totalDimension_;} ///< Get the (single) dimension for a unary system

  size_t getDimension(size_t i) const {return dimensions_[i];}

  void setDimensions(const Dimensions& dimensions) {dimensions_=dimensions; totalDimension_=product(dimensions);} ///< This will work only in the non-const case

private:
  Dimensions dimensions_;
  size_t totalDimension_ ;

};


/// dimensionality comparison for types derived from DimensionsBookkeeper \related DimensionsBookkeeper
template<int RANK>
inline bool
operator==(const DimensionsBookkeeper<RANK>& d1, const DimensionsBookkeeper<RANK>& d2)
{
  return blitz::all(d1.getDimensions()==d2.getDimensions());
}

/// dimensionality comparison for types derived from DimensionsBookkeeper \related DimensionsBookkeeper
template<int RANK>
inline bool
operator!=(const DimensionsBookkeeper<RANK>& d1, const DimensionsBookkeeper<RANK>& d2)
{
  return !(d1==d2);
}

#endif // CPPQEDCORE_QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED
