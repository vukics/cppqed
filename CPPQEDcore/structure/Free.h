// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_FREE_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_FREE_H_INCLUDED

#include "FreeFwd.h"

#include "DensityOperatorFwd.h"
#include "StateVectorFwd.h"

#ifndef FLAT_INCLUDE
#include "../quantumoperator/TridiagonalFwd.h"
#else //FLAT_INCLUDE
#include "TridiagonalFwd.h"
#endif //FLAT_INCLUDE
// Normally, structure is not allowed to depend on quantumoperator, here we make a small exception

#include "DynamicsBase.h"
#include "QuantumSystem.h"
#include "Types.h"

#include <boost/shared_ptr.hpp>

namespace structure {


/// Contains some typedefs for structures of arity 1 for convenience in defining free systems
/** \note The namespace was earlier called `free`, which in some situations created clashes with the global function of the same name declared in cstdlib */
namespace freesystem {

typedef quantumoperator::Tridiagonal<1> Tridiagonal; ///< unary Tridiagonal

typedef quantumdata::Types<1>::    StateVectorLow     StateVectorLow; ///< unary StateVectorLow
typedef quantumdata::Types<1>::DensityOperatorLow DensityOperatorLow; ///< unary DensityOperatorFwd

typedef quantumdata::LazyDensityOperator<1> LazyDensityOperator; ///< unary LazyDensityOperator

typedef quantumdata::    StateVector<1>     StateVector; ///< unary StateVector
typedef quantumdata::DensityOperator<1> DensityOperator; ///< unary DensityOperator

} // freesystem


/// In the language of the framework, a free system is a unary system (arity 1, `RANK=1`)
class Free : public QuantumSystem<1>, public DynamicsBase
{
public:
  typedef boost::shared_ptr<const Free> Ptr;

  /// A single dimension to initialise QuantumSystem`<1>` and the lists of real and complex name-frequency-multiplier tuples for DynamicsBase
  explicit Free(size_t dim, const RealFreqs& realFreqs=emptyRF, const ComplexFreqs& complexFreqs=emptyCF) : QuantumSystem<1>(dim), DynamicsBase(realFreqs,complexFreqs) {}

  Free(size_t dim, const ComplexFreqs& complexFreqs) : Free(dim,emptyRF,complexFreqs) {}
  Free(size_t dim, RealFreqsInitializer rf, ComplexFreqsInitializer cf={}) : Free(dim,RealFreqs(rf),ComplexFreqs(cf)) {}
  Free(size_t dim, ComplexFreqsInitializer cf) : Free(dim,{},cf) {}
  Free(size_t dim, RF rf, CF cf=CF()) : QuantumSystem<1>(dim), DynamicsBase(rf,cf) {}
  Free(size_t dim, CF cf) : QuantumSystem<1>(dim), DynamicsBase(cf) {}
  Free(size_t dim, RealFreqsInitializer rf, CF cf) : Free(dim,rf,{cf}) {}
  Free(size_t dim, RF rf, ComplexFreqsInitializer cf) : Free(dim,{rf},cf) {}

  /// \name Implementating inherited virtuals
  //@{
  /// Simply connects the pure virtual QuantumSystem::highestFrequency to the implementation DynamicsBase::highestFrequency.
  /**
   * An exception to the rule “never redefine an inherited non-virtual function” is taken because these functions are called the same in the two bases of Free,
   * which would otherwise create ambiguities.
   */
  double        highestFrequency (                ) const {return QuantumSystem<1>:: highestFrequency(  );}
  std::ostream& displayParameters(std::ostream& os) const {return QuantumSystem<1>::displayParameters(os);} ///< ”
  //@}

private:
  double         highestFrequency_v(                ) const final {return DynamicsBase::highestFrequency (  );}
  std::ostream& displayParameters_v(std::ostream& os) const final {return DynamicsBase::displayParameters(os);}

  std::ostream& displayMoreParameters(std::ostream& os) const final {return DynamicsBase::displayMoreParameters(os<<"# Dimension: "<<getDimension()<<std::endl);}

};


} // structure

#endif // CPPQEDCORE_STRUCTURE_FREE_H_INCLUDED
