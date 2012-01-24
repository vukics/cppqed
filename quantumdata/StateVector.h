// -*- C++ -*-
#ifndef _STATE_VECTOR_H
#define _STATE_VECTOR_H

#include "StateVectorFwd.h"

#include "DensityOperatorFwd.h"

#include "ArrayBase.h"
#include "DimensionsBookkeeper.h"
#include "LazyDensityOperator.h"
#include "Types.h"

#include "Operators.h"


namespace quantumdata {


template<int RANK1, int RANK2>
inline
const StateVector<RANK1+RANK2>
operator*(const StateVector<RANK1>&, const StateVector<RANK2>&);


template<int RANK>
const dcomp
braket(const StateVector<RANK>&, const StateVector<RANK>&);


// RANK is the number of quantum numbers


template<int RANK>
class StateVector 
  : public LazyDensityOperator<RANK>, 
    private ArrayBase<RANK>,
    private linalg::VectorSpace<StateVector<RANK> >
{
public:
  static const int N_RANK=RANK;

  typedef LazyDensityOperator<RANK> LDO_Base;
  typedef ArrayBase          <RANK>    ABase;

  typedef typename LDO_Base::Dimensions         Dimensions        ;

  typedef typename ABase::ArrayLow StateVectorLow;

  typedef typename Types<RANK>::DensityOperatorLow DensityOperatorLow;

  typedef typename LDO_Base::Idx Idx;

  using LDO_Base::getTotalDimension; using ABase::operator(); using ABase::vectorView;

  StateVector(const StateVectorLow& psi, ByReference) : LDO_Base(psi.shape()), ABase(psi) {}
  // By reference semantics! Simply furnishes an already existing StateVectorLow with a StateVector interface.
  // The trailing dummy argument is there only to make the user conscious of this fact.

  explicit StateVector(const Dimensions&, bool init=true);

  StateVector(const StateVector&); // use by value semantics (deep copy)

  // Default assignment is fine.

  // direct product
  template<int RANK2>
  StateVector(const StateVector<RANK2>&, const StateVector<RANK-RANK2>&);

  double   norm() const {return ABase::frobeniusNorm();}
  double renorm()                                     ;

  const DensityOperatorLow dyad(const StateVector&) const;
  const DensityOperatorLow dyad(                  ) const {return dyad(*this);}

  template<typename OTHER>
  StateVector& operator=(const OTHER& other) {operator()()=other; return *this;}
  // Together with the default assigment, this covers a lot of possibilities, including assignment from a StateVectorLow

  // naive operations for vector space

  StateVector& operator+=(const StateVector& psi) {ABase::operator+=(psi); return *this;}
  StateVector& operator-=(const StateVector& psi) {ABase::operator-=(psi); return *this;}

  const StateVector operator-() const {return -operator()();}
  const StateVector operator+() const {return *this;} // deep copy

  template<typename OTHER>
  StateVector& operator*=(const OTHER& dc) {ABase::operator*=(dc); return *this;}

  template<typename OTHER>
  StateVector& operator/=(const OTHER& dc) {ABase::operator/=(dc); return *this;}

  void addTo(DensityOperator<RANK>&) const;

protected:
  // virtual from LDO_Base
  const dcomp operator()(const Idx& i, const Idx& j) const {return operator()()(i)*conj(operator()()(j));}

};


template<int RANK1, int RANK2>
inline
const StateVector<RANK1+RANK2>
operator*(const StateVector<RANK1>& t1, const StateVector<RANK2>& t2)
{
  return StateVector<RANK1+RANK2>(t1,t2);
}


} // quantumdata

#include "impl/StateVector.tcc"

#endif // _STATE_VECTOR_H
