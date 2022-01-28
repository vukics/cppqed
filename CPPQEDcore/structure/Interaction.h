// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_INTERACTION_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_INTERACTION_H_INCLUDED

#include "DynamicsBase.h"
#include "Free.h"

#include "BlitzTiny.h"

#include <algorithm>

namespace structure {


/// Describes interaction of Free elements
/**
 * \tparam RANK Arity of the interaction. (Not greater than the arity of the full system Hilbert space.)
 * 
 * \note This is the simplest possible implementation allowing for interaction between frees only. If we are to achieve recursiveness,
 * i.e. the possibility of nesting composite systems into even more composite ones, we have to allow interaction between composites as well.
 * That is, between QuantumSystems in general. For this, Interaction should be even more templated taking compile-time vectors.
 * These specify between which quantum numbers of the subsystems the interaction acts. As many compile-time vectors are needed as the number of subsystems.
 * 
 * No inheritance from QuantumSystem because it does not make sense to simulate such an element as describes an interaction alone.
 * However, an interaction can have frequency-like parameters, hence the inheritance from DynamicsBase.
 * 
 */
template<int RANK>
class Interaction : public DynamicsBase, public DimensionsBookkeeper<RANK>
{
public:
  /// A tiny vector of shared pointers to the Free objects between which the interaction is defined
  /** \note The order of the Free objects is essential! (Cf. BinarySystem, Composite) */
  typedef std::array<FreePtr,RANK> Frees;

  typedef typename DimensionsBookkeeper<RANK>::Dimensions Dimensions;
  
  explicit Interaction(const Frees& frees,
                       const    RealFreqs&    realFreqs=emptyRF, 
                       const ComplexFreqs& complexFreqs=emptyCF)
    : DynamicsBase(realFreqs,complexFreqs), DimensionsBookkeeper<RANK>(extractDimensions(frees)), frees_(frees) {}

  auto operator[](size_t i) const {return frees_[i];}

private:
  static const Dimensions extractDimensions(const Frees& frees)
  {
    Dimensions res;
    std::transform(frees.begin(),frees.end(),res.begin(),[](auto f){return f->getTotalDimension();});
    return res;
  }
  
  const Frees frees_;

};


template<int RANK>
using InteractionPtr = std::shared_ptr<const Interaction<RANK>>;


} // structure

#endif // CPPQEDCORE_STRUCTURE_INTERACTION_H_INCLUDED
