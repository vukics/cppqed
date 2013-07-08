// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_INTERACTION_H_INCLUDED
#define STRUCTURE_INTERACTION_H_INCLUDED

#include "InteractionFwd.h"

#include "DynamicsBase.h"
#include "Free.h"

#include "BlitzTiny.h"

#include <boost/bind.hpp>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/copy.hpp>


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
 */
template<int RANK>
class Interaction : public DynamicsBase, public DimensionsBookkeeper<RANK>
{
public:
  typedef boost::shared_ptr<const Interaction> Ptr;

  /// A tiny vector of shared pointers to the Free objects between which the interaction is defined
  /** \note The order of the Free objects is essential! (Cf. BinarySystem, Composite) */
  typedef blitz::TinyVector<Free::Ptr,RANK> Frees;

  typedef typename DimensionsBookkeeper<RANK>::Dimensions Dimensions;
  
  explicit Interaction(const Frees& frees,
                       const    RealFreqs&    realFreqs=   RealFreqs(), 
                       const ComplexFreqs& complexFreqs=ComplexFreqs())
    : DynamicsBase(realFreqs,complexFreqs), DimensionsBookkeeper<RANK>(extractDimensions(frees)), frees_(frees) {}

    
  const Frees& getFrees() const {return frees_;}

private:
  static const Dimensions extractDimensions(const Frees& frees)
  {
    Dimensions res;
    boost::copy(frees | boost::adaptors::transformed(boost::bind(&Free::getTotalDimension,_1)), res.begin());
    return res;
  }
  
  const Frees frees_;

};


} // structure

#endif // STRUCTURE_INTERACTION_H_INCLUDED
