// -*- C++ -*-
#ifndef QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED
#define QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED

#include "DensityOperatorFwd.h"

#include "StateVectorFwd.h"

#include "ArrayBase.h"
#include "DimensionsBookkeeper.h"
#include "LazyDensityOperator.h"
#include "Types.h"

#include "Operators.h"


namespace quantumdata {


template<int RANK>
inline
const DensityOperator<RANK>
dyad(const StateVector<RANK>&, const StateVector<RANK>&);


template<int RANK1, int RANK2>
inline
const DensityOperator<RANK1+RANK2>
operator*(const DensityOperator<RANK1>&, const DensityOperator<RANK2>&);


template<int RANK>
inline
double frobeniusNorm(const DensityOperator<RANK>& rho) {return rho.frobeniusNorm();}


template<int RANK>
class DensityOperator
  : public LazyDensityOperator<RANK>,
    private ArrayBase<2*RANK>,
    private linalg::VectorSpace<DensityOperator<RANK> >
{
public:
  typedef LazyDensityOperator<  RANK> LDO_Base;
  typedef ArrayBase          <2*RANK>    ABase;

  typedef typename LDO_Base::Dimensions Dimensions;

  typedef typename LDO_Base::Idx Idx;

  typedef typename ABase::ArrayLow DensityOperatorLow;

  typedef linalg::CMatrix CMatrix;

  using ABase::operator(); using ABase::frobeniusNorm;

  /*
  DensityOperator() : Base() {}
  DensityOperatorBase() : LDO_Base(Dimensions()), ABase(DensityOperatorLow()), matrixView_() {}
  */

  DensityOperator(const DensityOperatorLow&, ByReference);

  explicit DensityOperator(const Dimensions&, bool init=true);

  explicit DensityOperator(const StateVector<RANK>&);

  DensityOperator(const DensityOperator&); // By value semantics (deep copy)

  // Default assignment doesn't work, because LazyDensityOperator is always purely constant (const DimensionsBookkeeper base)
  DensityOperator& operator=(const DensityOperator& rho) {ABase::operator=(rho()); return *this;}

  template<typename OTHER>
  DensityOperator& operator=(const OTHER& other) {operator()()=other; return *this;}
  
  // (multi-)matrix style indexing:
  const dcomp& operator()(const Idx& i, const Idx& j) const;
        dcomp& operator()(const Idx& i, const Idx& j) {return const_cast<dcomp&>(static_cast<const DensityOperator&>(*this)(i,j));}

  double   norm() const;
  double renorm()      ;

  const CMatrix matrixView() const {return blitzplusplus::binaryArray(operator()());}
  CMatrix       matrixView()       {return blitzplusplus::binaryArray(operator()());}

  // naive operations for vector space

  DensityOperator& operator+=(const DensityOperator& rho) {ABase::operator+=(rho); return *this;}
  DensityOperator& operator-=(const DensityOperator& rho) {ABase::operator-=(rho); return *this;}

  const DensityOperator operator-() const {return DensityOperator(-operator()());}
  const DensityOperator operator+() const {return *this;}

  template<typename OTHER>
  DensityOperator& operator*=(const OTHER& dc) {ABase::operator*=(dc); return *this;}

  template<typename OTHER>
  DensityOperator& operator/=(const OTHER& dc) {ABase::operator/=(dc); return *this;}

private:
  // virtual from LDO_Base
  const dcomp index(const Idx& i, const Idx& j) const {return operator()(i,j); }

};


template<int RANK>
void inflate(const TTD_DARRAY(1)&, DensityOperator<RANK>&, bool offDiagonals);


template<int RANK>
const DensityOperator<RANK>
densityOperatorize(const LazyDensityOperator<RANK>&);


template<int... SUBSYSTEM, int RANK>
const DensityOperator<mpl::size<tmptools::Vector<SUBSYSTEM...> >::value>
reduce(const LazyDensityOperator<RANK>&);


} // quantumdata


#endif // QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED