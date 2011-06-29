// -*- C++ -*-
#ifndef _FREE_H
#define _FREE_H

#include "FreeFwd.h"

#include "DensityOperatorFwd.h"
#include "StateVectorFwd.h"

#include "../quantumoperator/TridiagonalFwd.h"
// Normally, structure is not allowed to depend on quantumoperator, here we make a small exception

#include "DynamicsBase.h"
#include "QuantumSystem.h"
#include "Types.h"

namespace structure {


namespace free {

typedef quantumoperator::Tridiagonal<1> Tridiagonal;

typedef quantumdata::Types<1>::    StateVectorLow     StateVectorLow;
typedef quantumdata::Types<1>::DensityOperatorLow DensityOperatorLow;

typedef quantumdata::LazyDensityOperator<1> LazyDensityOperator;

typedef quantumdata::    StateVector<1>     StateVector;
typedef quantumdata::DensityOperator<1> DensityOperator;

} // free


class Free : public QuantumSystem<1>, public DynamicsBase
{
public:
  explicit Free(size_t, const RealFreqs& =RealFreqs(), const ComplexFreqs& =ComplexFreqs());

  double highestFrequency (                ) const {return DynamicsBase::highestFrequency (  );}
  void   displayParameters(std::ostream& os) const {return DynamicsBase::displayParameters(os);}

private:
  void displayMoreParameters(std::ostream&) const;

};


} // structure

#endif // _FREE_H
