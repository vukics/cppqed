// -*- C++ -*-
#ifndef STRUCTURE_FREE_H_INCLUDED
#define STRUCTURE_FREE_H_INCLUDED

#include "FreeFwd.h"

#include "DensityOperatorFwd.h"
#include "StateVectorFwd.h"

#include "../quantumoperator/TridiagonalFwd.h"
// Normally, structure is not allowed to depend on quantumoperator, here we make a small exception

#include "DynamicsBase.h"
#include "QuantumSystem.h"
#include "Types.h"

#include <boost/shared_ptr.hpp>

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
  typedef boost::shared_ptr<const Free> SmartPtr;

  explicit Free(size_t, const RealFreqs& =RealFreqs(), const ComplexFreqs& =ComplexFreqs());

  double highestFrequency (                ) const {return DynamicsBase::highestFrequency (  );}
  void   displayParameters(std::ostream& os) const {return DynamicsBase::displayParameters(os);}

private:
  void displayMoreParameters(std::ostream&) const;

};


} // structure

#endif // STRUCTURE_FREE_H_INCLUDED
