// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_FREE_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_FREE_H_INCLUDED

#include "DensityOperator.h"
#include "DynamicsBase.h"
#include "QuantumSystem.h"
#include "Types.h"


namespace structure {


/// Contains some typedefs for structures of arity 1 for convenience in defining free systems
/** \note The namespace was earlier called `free`, which in some situations created clashes with the global function of the same name declared in cstdlib */
namespace freesystem {

typedef quantumdata::StateVectorLow<1> StateVectorLow; ///< unary StateVectorLow

typedef quantumdata::DensityOperatorLow<1> DensityOperatorLow; ///< unary DensityOperatorFwd

typedef quantumdata::LazyDensityOperator<1> LazyDensityOperator; ///< unary LazyDensityOperator

typedef quantumdata::    StateVector<1>     StateVector; ///< unary StateVector
typedef quantumdata::DensityOperator<1> DensityOperator; ///< unary DensityOperator

using ::structure::Averages; using ::structure::Rates;

} // freesystem


/// In the language of the framework, a free system is a unary system (arity 1, `RANK=1`)
class Free : public QuantumSystem<1>, public DynamicsBase
{
public:
  /// A single dimension to initialise QuantumSystem`<1>` and the lists of real and complex name-frequency-multiplier tuples for DynamicsBase
  explicit Free(size_t dim, const RealFreqs& realFreqs={}, const ComplexFreqs& complexFreqs={}) : QuantumSystem<1>(dim), DynamicsBase(realFreqs,complexFreqs) {}

  /// \name Implementating inherited virtuals
  //@{
  /// Simply connects the pure virtual QuantumSystem::highestFrequency to the implementation DynamicsBase::highestFrequency.
  /**
   * An exception to the rule “never redefine an inherited non-virtual function” is taken because these functions are called the same in the two bases of Free,
   * which would otherwise create ambiguities.
   */
  using DynamicsBase::highestFrequency;
  using DynamicsBase::streamParameters;
  //@}

private:
  double         highestFrequency_v(                ) const final {return DynamicsBase::highestFrequency (  );}
  std::ostream& streamParameters_v(std::ostream& os) const final {return DynamicsBase::streamParameters(os);}

  std::ostream& streamMoreParameters(std::ostream& os) const final {return DynamicsBase::streamMoreParameters(os<<"Dimension: "<<getDimension()<<std::endl);}

};


using FreePtr = std::shared_ptr<const Free>;


} // structure

#endif // CPPQEDCORE_STRUCTURE_FREE_H_INCLUDED
