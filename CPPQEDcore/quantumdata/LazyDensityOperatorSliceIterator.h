// -*- C++ -*-
/// \briefFile{Defines quantumdata::ldo::DiagonalIterator}
#ifndef QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_H_INCLUDED

#include "LazyDensityOperatorFwd.h"

#include <boost/operators.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/mpl/size.hpp>


namespace mpl=boost::mpl;


namespace quantumdata {

/// Comprises tools related to <strong>L</strong>azy<strong>D</strong>ensity<strong>O</strong>perator, like DiagonalIterator
namespace ldo {
  
/** \page slicinganldo Slicing a LazyDensityOperator
 * 
 * \tableofcontents
 * 
 * Analogously to \ref multiarrayconcept "slicing state vectors", it is also necessary to slice LazyDensityOperator objects
 * because for calculating quantum expectation values of subsystem-observables (e.g. in Composite objects),
 * the partial-trace density operator is needed.
 * 
 * For the partial trace, however, only such elements of the full density operator are needed as are diagonal in the indices *not* belonging
 * to the given subsystem (dummy indices). This is the reason why the tool performing the iteration is called DiagonalIterator.
 *
 * Slicing is fully recursive, a sliced LazyDensityOperator (usually obtained by dereferencing a DiagonalIterator) can be further sliced.
 * 
 * Notes on implementation {#slicinganldoimplementation}
 * =======================
 * 
 * Difficulty: LazyDensityOperator is an abstract interface, and the organization of its data (and hence the actual procedure of slicing)
 * varies along implementations.
 * 
 * \image html ldoDiagonalIterator.png
 * 
 * Implemented using a classical inheritance-based strategy idiom, together with both compile-time and run-time implementation selection
 * (similarly to blitzplusplus::basi::Iterator, a special implementation (DiagonalIterator::DI_ImplSpecial) is needed when the size
 * of the compile-time vector `V` equals `RANK`).
 * 
 * The run-time polymorphy of LazyDensityOperator necessitates run-time implementation selection.
 * This happens when the actual implementation (`DI_SV_Impl` or `DI_DO_Impl`, or their non-orthogonal counterparts) is chosen 
 * depending on whether the given LazyDensityOperator is a StateVector or a DensityOperator (or their non-orthogonal counterparts).
 * 
 */

#define INPUT_IteratorHelper boost::input_iterator_helper<DiagonalIterator<RANK,V>,const LazyDensityOperator<mpl::size<V>::value> >


/// Iterator for slices of a LazyDensityOperator that are diagonal in the dummy indices
/**
 * Cf. \ref slicinganldo "rationale"
 * 
 * \tparam RANK arity of the full (unsliced) Hilbert space
 * \tparamV
 * 
 * It's inherently a const iterator since LazyDensityOperator is immutable.
 * Models an [InputIterator](http://www.cplusplus.com/reference/std/iterator/InputIterator/), implemented with the help of
 * \refBoostConstruct{input_iterator_helper,utility/operators.htm#iterator} from Boost.Operator.
 * 
 */
template<int RANK, typename V>
class DiagonalIterator 
  : public INPUT_IteratorHelper 
{
private:
  typedef INPUT_IteratorHelper Base;

#undef INPUT_IteratorHelper

  typedef typename Base::value_type LazyDensityOperatorRes; ///< Equivalent to `const LazyDensityOperator<mpl::size<V>::value>`

public:
  /// Constructor
  /**
   * Similarly to blitzplusplus::basi::Iterator, it can be initialised either to the beginning or to the end of the sequence.
   * \tparam IS_END governs the end-ness
   */
  template<bool IS_END>
  DiagonalIterator(const LazyDensityOperator<RANK>& ldo, mpl::bool_<IS_END>);
  
  /// \name Necessary members of an input iterator
  //@{
  DiagonalIterator& operator++() {impl_->doIncrement(); return *this;} ///< Immediately delegated to the implementation

  const LazyDensityOperatorRes& operator*() const {return impl_->dereference();} ///< ”
  
  bool operator==(const DiagonalIterator& other) const {return impl_->isEqualTo(*other.impl_);} ///< ”
  //@}
  
  static const bool IS_SPECIAL=(RANK==mpl::size<V>::value); ///< Signifies whether the special implementation is needed

  class DI_ImplSpecial; ///< Special implementation in the case `RANK=mpl::size<V>::value`
  class DI_Impl       ; ///< Base class for non-special implementations

  /// Pointer to implementation
  /**
   * We are using the classical inheritance-based pointer-to-implementation technique, together with some compile-time dispatching.
   * \refBoostConstruct{eval_if,mpl/doc/refmanual/eval-if.html} from Boost.MPL here guarantees that `DI_ImplSpecial` gets instantiated
   * only in the special case
   */
  typedef boost::shared_ptr<typename mpl::eval_if_c<IS_SPECIAL,
                                                    mpl::identity<DI_ImplSpecial>,
                                                    mpl::identity<DI_Impl       > 
                                                    >::type
                           > Impl;

private:
  Impl impl_;
  
};

} // ldo

} // quantumdata


#endif // QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_H_INCLUDED
