// -*- C++ -*-
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_INTERACTION_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_INTERACTION_H_INCLUDED

#include "InteractionFwd.h"

#include "DynamicsBase.h"
#include "Free.h"

#include "BlitzTiny.h"
#include "SmartPtr.h"

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
 * No inheritance from QuantumSystem because it does not make sense to simulate such an element as describes an interaction alone.
 * However, an interaction can have frequency-like parameters, hence the inheritance from DynamicsBase.
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

  template<typename F> using FreesTemplate=blitz::TinyVector<F,RANK>; ///< `F` must be convertible to Free::Ptr by cpputils::sharedPointerize
  
  typedef typename DimensionsBookkeeper<RANK>::Dimensions Dimensions;
  
  /// \todo Why not a variadic constructor here, taking a pack of Free instances?
  explicit Interaction(const Frees& frees,
                       const    RealFreqs&    realFreqs=emptyRF, 
                       const ComplexFreqs& complexFreqs=emptyCF)
    : DynamicsBase(realFreqs,complexFreqs), DimensionsBookkeeper<RANK>(extractDimensions(frees)), frees_(frees) {}

  Interaction(const Frees& frees, const ComplexFreqs& complexFreqs) : Interaction(frees,emptyRF,complexFreqs) {}
  Interaction(const Frees& frees, RealFreqsInitializer rf, ComplexFreqsInitializer cf={}) : Interaction(frees,RealFreqs(rf),ComplexFreqs(cf)) {}
  Interaction(const Frees& frees, ComplexFreqsInitializer cf) : Interaction(frees,RealFreqsInitializer(),cf) {}
  Interaction(const Frees& frees, RF rf, CF cf=CF()) : Interaction(frees,RealFreqsInitializer{rf}, cf==CF() ? ComplexFreqsInitializer{} : ComplexFreqsInitializer{cf}) {}
  Interaction(const Frees& frees, CF cf) : Interaction(frees,ComplexFreqsInitializer{cf}) {}
  Interaction(const Frees& frees, RealFreqsInitializer rf, CF cf) : Interaction(frees,rf,{cf}) {}
  Interaction(const Frees& frees, RF rf, ComplexFreqsInitializer cf) : Interaction(frees,{rf},cf) {}

  template<typename F0, typename F1, typename... F_REST, typename RF_type, typename CF_type>
  Interaction(F0&& f0, F1&& f1, F_REST&&... f_rest, const RealFreqs& realFreqs=emptyRF, const ComplexFreqs& complexFreqs=emptyCF)
    :  Interaction(Frees(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f_rest)...),realFreqs,complexFreqs) {}

  template<typename F0, typename F1, typename... F_REST>
  Interaction(F0&& f0, F1&& f1, F_REST&&... f_rest, const ComplexFreqs& complexFreqs) : Interaction(Frees(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f_rest)...),
                                                                                                    complexFreqs) {}

  template<typename F0, typename F1, typename... F_REST>
  Interaction(F0&& f0, F1&& f1, F_REST&&... f_rest, RealFreqsInitializer rf, ComplexFreqsInitializer cf={}) : Interaction(Frees(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f_rest)...),
                                                                                                                          rf,cf) {}

  template<typename F0, typename F1, typename... F_REST>
  Interaction(F0&& f0, F1&& f1, F_REST&&... f_rest, ComplexFreqsInitializer cf) : Interaction(Frees(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f_rest)...),
                                                                                              cf) {}

  template<typename F0, typename F1, typename... F_REST>
  Interaction(F0&& f0, F1&& f1, F_REST&&... f_rest, RF rf, CF cf=CF()) : Interaction(Frees(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f_rest)...),
                                                                                     rf,cf) {}

  template<typename F0, typename F1, typename... F_REST>
  Interaction(F0&& f0, F1&& f1, F_REST&&... f_rest, CF cf) : Interaction(Frees(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f_rest)...),
                                                                         cf) {}

  template<typename F0, typename F1, typename... F_REST>
  Interaction(F0&& f0, F1&& f1, F_REST&&... f_rest, RealFreqsInitializer rf, CF cf) : Interaction(Frees(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f_rest)...),
                                                                                                  rf,cf) {}

  template<typename F0, typename F1, typename... F_REST>
  Interaction(F0&& f0, F1&& f1, F_REST&&... f_rest, RF rf, ComplexFreqsInitializer cf) : Interaction(Frees(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f_rest)...),
                                                                                                     rf,cf) {}

  const Frees& getFrees() const {return frees_;}

protected:
  /// Shared-pointerizes the elements passed as Frees. \todo How to shared-pointerize an rvalue reference?
  class FreesProxy
  {
  public:
    template<typename... Pack>
    FreesProxy(Pack&&... pack) : frees_(cpputils::sharedPointerize(pack)...) {}
    
    operator const Frees&() const {return frees_;}
    
  private:
    const Frees frees_;
  };

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

#endif // CPPQEDCORE_STRUCTURE_INTERACTION_H_INCLUDED
