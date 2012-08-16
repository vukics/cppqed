// -*- C++ -*-
#ifndef QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED
#define QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED

#include "DimensionsBookkeeperFwd.h"

#include "Exception.h"
#include "TMP_Tools.h"
#include "BlitzTiny.h"

namespace mpl=boost::mpl;


struct DimensionalityMismatchException : cpputils::Exception {};


template<int RANK, bool IS_CONST>
class DimensionsBookkeeper
{
public:
  static const int N_RANK=RANK;

  typedef TTD_EXTTINY(RANK) Dimensions;

  explicit DimensionsBookkeeper(mpl::bool_<IS_CONST> =mpl::false_()) : dimensions_(), totalDimension_() {}

  explicit DimensionsBookkeeper(const Dimensions& dimensions) : dimensions_(dimensions), totalDimension_(product(dimensions)) {}

  const Dimensions& getDimensions    () const {return     dimensions_;}
  size_t            getTotalDimension() const {return totalDimension_;}

  size_t getDimension(mpl::int_<RANK> =mpl::int_<1>()) const {return totalDimension_;} 

  void setDimensions(const Dimensions& dimensions) {dimensions_=dimensions; totalDimension_=product(dimensions);}

private:
  typename tmptools::ConditionalAddConst<Dimensions,IS_CONST>::type      dimensions_;
  typename tmptools::ConditionalAddConst<size_t    ,IS_CONST>::type totalDimension_ ;

};


template<int RANK, bool IS_CONST1, bool IS_CONST2>
inline bool
operator==(const DimensionsBookkeeper<RANK,IS_CONST1>& d1, const DimensionsBookkeeper<RANK,IS_CONST2>& d2)
{
  return blitz::all(d1.getDimensions()==d2.getDimensions());
}


template<int RANK, bool IS_CONST1, bool IS_CONST2>
inline bool
operator!=(const DimensionsBookkeeper<RANK,IS_CONST1>& d1, const DimensionsBookkeeper<RANK,IS_CONST2>& d2)
{
  return !(d1==d2);
}


#endif // QUANTUMDATA_DIMENSIONSBOOKKEEPER_H_INCLUDED
