// -*- C++ -*-
/// \briefFileDefault
#ifndef QUANTUMDATA_STATEVECTOR_H_INCLUDED
#define QUANTUMDATA_STATEVECTOR_H_INCLUDED

#include "StateVectorFwd.h"

#include "DensityOperatorFwd.h"

#include "ArrayBase.h"
#include "DimensionsBookkeeper.h"
#include "LazyDensityOperator.h"
#include "Types.h"

#include "Operators.h"


namespace quantumdata {

/** \page quantumdatahighlevel High-level data structures
 * 
 * The StateVector and DensityOperator classes (and their non-orthogonal counterparts) exist for two main reasons:
 * - As interfaces to Types::StateVectorLow and Types::DensityOperatorLow, respectively, which are more convenient to use on higher levels of the framework,
 * e.g. in scripts and in quantum trajectories. This is especially because while `blitz::Array` uses by-reference copy semantics and expression templates,
 * these classes use the more usual by-value copy semantics and normal semantics for arithmetic operations. This means, however, 
 * that copying and arithmetics should be used judiciously, if possible only in the startup phase of simulations.
 * - As implementations of the LazyDensityOperator interface.
 * 
 */

/// Creates the direct product, relying on the direct-product constructor
template<int RANK1, int RANK2>
inline
const StateVector<RANK1+RANK2>
operator*(const StateVector<RANK1>&, const StateVector<RANK2>&);


/// Calculates the inner product, relying on StateVector::vectorView
template<int RANK>
const dcomp
braket(const StateVector<RANK>&, const StateVector<RANK>&);


/// State vector of arbitrary arity
/**
 * Cf. \ref quantumdatahighlevel "rationale"
 * 
 * \tparamRANK
 * 
 * The inheritance of StateVector from linalg::VectorSpace provides for a lot of free-standing helpers describing vector-space algebra.
 * These are all naively based on the arithmetic member functions like StateVector::operator+=, StateVector::operator*=, etc.
 * 
 */
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

  /// \name Construction, assignment
  //@{
    /// Constructs the class in such a way that the underlying data reference the same data as `psi`.
    /**
    * Simply furnishes an already existing StateVectorLow with a StateVector interface.
    * 
    * \note By-reference semantics! Since everywhere else the class represents by-value semantics, some care is needed with this constructor’s use.
    * For this reason, the tagging dummy class ByReference is introduced, to make the user conscious of what the semantics is.
    */
  StateVector(const StateVectorLow& psi, ByReference) : LDO_Base(psi.shape()), ABase(psi) {}

  explicit StateVector(const Dimensions&, bool init=true); ///< Constructs the class with a newly allocated chunk of memory, which is initialized only if `init` is `true`.

  StateVector(const StateVector&); ///< Copy constructor using by value semantics, that is, deep copy.

    /// Constructs the class as the direct product of `psi1` and `psi2`, whose arities add up to `RANK`.
    /**
    * The implementation relies on blitzplusplus::concatenateTinies and blitzplusplus::doDirect.
    * \tparam RANK2 the arity of one of the operands
    */
  template<int RANK2>
  StateVector(const StateVector<RANK2>&, const StateVector<RANK-RANK2>&);
  
    /// Assignment with by-value semantics.
    /** Default assignment doesn't work, because LazyDensityOperator is always purely constant (const DimensionsBookkeeper base). */
  StateVector& operator=(const StateVector& sv) {ABase::operator=(sv()); return *this;}

    /// Mixed-mode assignment with by-value semantics
    /**
    * The standard assignment and the templated assignment together cover a lot of possibilities, including also assignment from a StateVectorLow,
    * but for example also from a DArray<RANK>, or just a const c-number. (Can be assigned from anything a CArray<RANK> can be assigned from.)
    * 
    * \tparam OTHER the “other” type in mixed mode
    */
  template<typename OTHER>
  StateVector& operator=(const OTHER& other) {operator()()=other; return *this;}
  //@}
  
  /// \name Metric
  //@{
    /// Both functions return the norm \f$\norm\Psi\f$, but the latter one also renormalises
    /** Implemented in terms of ArrayBase::frobeniusNorm. */
  double   norm() const {return ABase::frobeniusNorm();}
  double renorm()                                     ;
  //@}
  
  /// \name Dyad
  //@{
    /// Both functions form a dyad, the second one with the object itself
    /** This is a rather expensive operation, implemented in terms of blitzplusplus::doDirect. */
  const DensityOperatorLow dyad(const StateVector&) const;
  const DensityOperatorLow dyad(                  ) const {return dyad(*this);}
  //@}

  /// Adds a dyad of the present object to `densityOperator`
  /**
   * This is done without actually forming the dyad in memory (so that this is not implemented in terms of StateVector::dyad).
   * This is important in situations when an average density operator is needed from an ensemble of state vectors, an example being quantumtrajectory::EnsembleMCWF.
   */
  void addTo(DensityOperator<RANK>& densityOperator) const;

  /// \name Naive vector-space operations
  //@{
  StateVector& operator+=(const StateVector& psi) {ABase::operator+=(psi); return *this;}
  StateVector& operator-=(const StateVector& psi) {ABase::operator-=(psi); return *this;}

  const StateVector operator-() const {StateVector res(this->getDimensions(),false); res()=-this->operator()(); return res;}
  const StateVector operator+() const {return *this;} // deep copy
  //@}

  /// \name Naive vector-space operations allowing also for mixed-mode arithmetics
  //@{
    /** \tparam OTHER the “other” type in mixed mode */
  template<typename OTHER>
  StateVector& operator*=(const OTHER& dc) {ABase::operator*=(dc); return *this;}

  template<typename OTHER>
  StateVector& operator/=(const OTHER& dc) {ABase::operator/=(dc); return *this;}
  //@}
  
private:
  const dcomp index(const Idx& i, const Idx& j) const {return operator()()(i)*conj(operator()()(j));} ///< This function implements the LazyDensityOperator interface in a dyadic-product way.

};


template<int RANK1, int RANK2>
inline
const StateVector<RANK1+RANK2>
operator*(const StateVector<RANK1>& t1, const StateVector<RANK2>& t2)
{
  return StateVector<RANK1+RANK2>(t1,t2);
}


} // quantumdata

#endif // QUANTUMDATA_STATEVECTOR_H_INCLUDED
