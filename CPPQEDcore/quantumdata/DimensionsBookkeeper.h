// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED

#include "DimensionsBookkeeperFwd.h"

#include "Exception.h"
#include "TMP_Tools.h"
#include "BlitzTiny.h"

namespace mpl=boost::mpl;


/// Thrown in the case of dimensionality mismatch of constructs over the same Hilbert space
class DimensionalityMismatchException : cpputils::Exception {};


/// Stores and manipulates dimensions of constructs over composite Hilbert spaces of arbitrary arity
/**
 * It can be either constant or non-constant depending on the second template parameter.
 * 
 * \tparamRANK
 * \tparam IS_CONST  Governs the constness
 */
template<int RANK, bool IS_CONST>
class DimensionsBookkeeper
{
public:
  static const int                    N_RANK=RANK; ///< Arity of the Hilbert space
  static const int DIMESIONS_BOOKKEEPER_RANK=RANK; ///< Ditto (to break ambiguity if a class is derived from another base featuring `N_RANK`).

  typedef ExtTiny<RANK> Dimensions; ///< The dimensions as a static vector of size N_RANK

  /// Constructor usable only in the `IS_CONST=false` case
  /**
   * The aim of the dummy argument with a default value – which creates a nonsensical function signature for `IS_CONST=true` – is that 
   * this constructor only compile for `IS_CONST=false` because it is only in the non-constant case that we allow default construction of the class.
   * 
   * Since from a template only such parts are compiled as are actually used, a client can use the class in the case 
   * `IS_CONST=true` without problems, getting a compile-time error only when trying to default-construct such an object.
   */
  explicit DimensionsBookkeeper(mpl::bool_<IS_CONST> =mpl::false_())
    : dimensions_(), totalDimension_() {}

  explicit DimensionsBookkeeper(const Dimensions& dimensions)
    : dimensions_(dimensions), totalDimension_(product(dimensions)) {} ///< Standard constructor usable also in the `IS_CONST=true` case

  const Dimensions& getDimensions    () const {return     dimensions_;} ///< Get the Dimensions vector
  size_t            getTotalDimension() const {return totalDimension_;} ///< Get the total dimension of a system of arbitrary arity

  size_t getDimension(mpl::int_<RANK> =mpl::int_<1>()) const {return totalDimension_;} ///< Get the (single) dimension for a unary system

  size_t getDimension(size_t i) const {return dimensions_[i];}

  void setDimensions(const Dimensions& dimensions) {dimensions_=dimensions; totalDimension_=product(dimensions);} ///< This will work only in the non-const case

private:
  typename tmptools::ConditionalAddConst<Dimensions,IS_CONST>::type      dimensions_;
  typename tmptools::ConditionalAddConst<size_t    ,IS_CONST>::type totalDimension_ ;

};


/// dimensionality comparison for types derived from DimensionsBookkeeper \related DimensionsBookkeeper
template<int RANK, bool IS_CONST1, bool IS_CONST2>
inline bool
operator==(const DimensionsBookkeeper<RANK,IS_CONST1>& d1, const DimensionsBookkeeper<RANK,IS_CONST2>& d2)
{
  return blitz::all(d1.getDimensions()==d2.getDimensions());
}

/// dimensionality comparison for types derived from DimensionsBookkeeper \related DimensionsBookkeeper
template<int RANK, bool IS_CONST1, bool IS_CONST2>
inline bool
operator!=(const DimensionsBookkeeper<RANK,IS_CONST1>& d1, const DimensionsBookkeeper<RANK,IS_CONST2>& d2)
{
  return !(d1==d2);
}

#endif // CPPQEDCORE_QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED
