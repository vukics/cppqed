// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_TCC_INCLUDED

#include "DensityOperator.h"

#include "BlitzArrayExtensions.tcc"
#include "MultiIndexIterator.tcc"


namespace quantumdata {


template<int RANK>
inline
const DensityOperator<RANK>
dyad(const StateVector<RANK>& sv1, const StateVector<RANK>& sv2)
{
  return DensityOperator<RANK>(sv1.dyad(sv2),byReference);
}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const StateVector<RANK>& psi) 
  : LDO_Base(psi.getDimensions()),
    ABase(psi.dyad())
{}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const DensityOperatorLow& rho, ByReference) 
  : LDO_Base(blitzplusplus::halfCutTiny(rho.shape())),
    ABase(rho)
{}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const Dimensions& dimensions, bool init)
  : LDO_Base(dimensions),
    ABase(DensityOperatorLow(blitzplusplus::concatenateTinies(dimensions,dimensions)))
{
  if (init) *this=0;
}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const DensityOperator& rho) 
  : LDO_Base(rho.getDimensions()),
    ABase(rho.getArray().copy())
{}


template<int RANK>
DensityOperator<RANK>::DensityOperator(DensityOperator&& rho) 
  : LDO_Base(rho.getDimensions()),
    ABase(rho.getArray())
{}


template<int RANK>
double
DensityOperator<RANK>::norm() const
{
  using blitz::tensor::i;
  const linalg::CMatrix m(matrixView());
  return real(sum(m(i,i)));
}


template<int RANK>
double
DensityOperator<RANK>::renorm()
{
  double trace=norm();
  *this/=trace;
  return trace;
}



template<int RANK>
const dcomp&
DensityOperator<RANK>::indexWithTiny(const Idx& i, const Idx& j) const
{
  return getArray()(blitzplusplus::concatenateTinies<int,int,RANK,RANK>(i,j));
}



template<int RANK>
void inflate(const DArray<1>& flattened, DensityOperator<RANK>& rho, bool offDiagonals)
{
  using mathutils::sqr;

  const size_t dim=rho.getTotalDimension();
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(rho.getDimensions()-1,cpputils::mii::begin);
  
  size_t idx=0;

  // Diagonal
  for (Iterator i(etalon); idx<dim; ++i)
    rho(*i)(*i)=flattened(idx++);
  
  // OffDiagonal
  if (offDiagonals)
    for (Iterator i(etalon); idx<mathutils::sqr(dim); ++i)
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j, idx+=2) {
        dcomp matrixElement(rho(*i)(*j)=dcomp(flattened(idx),flattened(idx+1)));
        rho(*j)(*i)=conj(matrixElement);
      }

}


template<int RANK>
const DensityOperator<RANK>
densityOperatorize(const LazyDensityOperator<RANK>& matrix)
{
  DensityOperator<RANK> res(matrix.getDimension());
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(matrix.getDimensions()-1,cpputils::mii::begin);
  
  for (Iterator i(etalon); i!=etalon.getEnd(); ++i) {
    res(*i)(*i)=matrix(*i);
    for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j) {
      dcomp matrixElement(res(*i)(*j)=matrix(*i)(*j));
      res(*j)(*i)=conj(matrixElement);
    }
  }

  return res;
  
}


template<int... SUBSYSTEM, int RANK>
const DensityOperator<mpl::size<tmptools::Vector<SUBSYSTEM...> >::value>
reduce(const LazyDensityOperator<RANK>& matrix)
{
  static const int RES_ARITY=mpl::size<tmptools::Vector<SUBSYSTEM...> >::value;
  return partialTrace<tmptools::Vector<SUBSYSTEM...>,DensityOperator<RES_ARITY> >(matrix,densityOperatorize<RES_ARITY>);
}


} // quantumdata

#endif // CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_TCC_INCLUDED
