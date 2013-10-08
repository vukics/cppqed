// -*- C++ -*-
/// \briefFileDefault
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


/// Contains some typedefs for structures of arity 1 for convenience in defining free systems
namespace free {

typedef quantumoperator::Tridiagonal<1> Tridiagonal;

typedef quantumdata::Types<1>::    StateVectorLow     StateVectorLow;
typedef quantumdata::Types<1>::DensityOperatorLow DensityOperatorLow;

typedef quantumdata::LazyDensityOperator<1> LazyDensityOperator;

typedef quantumdata::    StateVector<1>     StateVector;
typedef quantumdata::DensityOperator<1> DensityOperator;

} // free


/// In the language of the framework, a free system is a unary system (arity 1, `RANK=1`)
class Free : public QuantumSystem<1>, public DynamicsBase
{
public:
  typedef boost::shared_ptr<const Free> Ptr;

  explicit Free(size_t, const RealFreqs& =RealFreqs(), const ComplexFreqs& =ComplexFreqs()); ///< A single dimension to initialise QuantumSystem`<1>` and the lists of real and complex name-frequency-multiplier tuples for DynamicsBase

  /// \name Implementating inherited virtuals
  /// Simply connect the pure virtual QuantumSystem::highestFrequency to the implementation DynamicsBase::highestFrequency.
  /**
   * An exception to the rule “never redefine an inherited non-virtual function” is taken because these functions are called the same in the two bases of Free,
   * which would otherwise create ambiguities.
   */
  // @{
  double        highestFrequency (                ) const {return  highestFrequency_v(  );}
  std::ostream& displayParameters(std::ostream& os) const {return displayParameters_v(os);}
  // @}

private:
  double         highestFrequency_v(                ) const {return DynamicsBase::highestFrequency (  );}
  std::ostream& displayParameters_v(std::ostream& os) const {return DynamicsBase::displayParameters(os);}

  std::ostream& displayMoreParameters(std::ostream&) const;

};


} // structure

#endif // STRUCTURE_FREE_H_INCLUDED
