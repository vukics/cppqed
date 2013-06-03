// -*- C++ -*-
/// \briefFileDefault
#ifndef QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED

#include "LazyDensityOperatorFwd.h"

#include "DimensionsBookkeeper.h"

#include "BlitzArray.h"
#include "ComplexExtensions.h"

#include <boost/shared_ptr.hpp>

#include <boost/mpl/if.hpp>


namespace mpl=boost::mpl;


/// Comprises classes representing the state of composite quantum systems and providing various interfaces to manipulate this data
/** Some of its most important classes fit into a single class rooted in the virtual interface LazyDensityOperator. */
namespace quantumdata {


/// The primary tool for performing \ref slicinganldo "slice iterations"
/**
 * On higher levels of the framework (cf. eg. BinarySystem, Composite), this function is used exclusively for performing LazyDensityOperator slice iteration.
 * 
 * \tparamV
 * \tparam T an arithmetic type which must be default-constructible
 * \tparam F a callable type with signature `const T(const typename ldo::DiagonalIterator<RANK,V>::value_type&)`. This signature is equivalent to `const T(const LazyDensityOperator<mpl::size<V>::value>&)`.
 * 
 */
template<typename V, typename T, int RANK, typename F>
const T
partialTrace(const LazyDensityOperator<RANK>&, F function);


/// Common interface for calculating quantum averages
/**
 * In a quantum-simulation framework, users should be able to write code for calculating quantum expectation values from quantumdata, 
 * independently of whether this data is represented by state vectors or density operators, in orthogonal or non-orthogonal bases. 
 * One obvious solution is relying on the formula
 * \f[\avr{A}=\Tr{A\rho}\f]
 * (\f$A\f$ being an observable and \f$\rho\f$ the density operator of the system), to write code only for the density-operator case,
 * and fall back to this in the state-vector case as well, by calculating a dyad from the state vector. This is, however, extremely wasteful, 
 * since usually not all the matrix elements of \f$\rho\f$ are needed for calculating the average, furthermore, for large dimensionality,
 * this solution may become outright unaffordable in terms of memory: for large systems, we may afford to store \f$\ket\Psi\f$, but not \f$\ket\Psi\bra\Psi\f$.
 * 
 * The solution adopted for this problem in the framework is represented by LazyDensityOperator, which provides a common interface 
 * for all the four cases StateVector, DensityOperator, and their non-orthogonal counterparts, to calculate quantum averages from their data.
 * The “laziness” amounts to that in the case of state vectors, only those elements of the density operator are calculated that are actually asked for.
 * 
 * This class is totally immutable, all its member functions being constant.
 * 
 * \tparamRANK
 */
template<int RANK> 
class LazyDensityOperator 
  : public DimensionsBookkeeper<RANK,true>
{
public:
  typedef boost::shared_ptr<const LazyDensityOperator> Ptr; ///< Many class templates in the framework define shared pointers to their own types, in a template-metafunction like manner
  
  typedef DimensionsBookkeeper<RANK,true> Base;

  typedef typename Base::Dimensions Dimensions; ///< Inherited from DimensionsBookkeeper

  /// The type used for indexing the “rows” and the “columns”.
  /** Just an integer (index) if `RANK=1`, otherwise a tiny vector of integers (multi-index). */
  typedef typename mpl::if_c<(RANK==1),int,TTD_IDXTINY(RANK)>::type Idx; 

  virtual ~LazyDensityOperator() {}

  const dcomp operator()(const Idx& i, const Idx& j) const {return index(i,j);} ///< The indexing function calling a purely virtual function according to the non-virtual interface idiom

  double operator()(const Idx& i) const {return real((*this)(i,i));} ///< An inline function for conveniently addressing the diagonal elements

  /// \name Slicing-related functionality
  //@{
    /// Return the ldo::DiagonalIterator corresponding to the beginning/end of the sequence of slices defined by `V` 
    /** Cf. \ref slicinganldo "rationale"
     * \tparamV
     */
  template<typename V>
  const ldo::DiagonalIterator<RANK,V> begin() const; 

  template<typename V>
  const ldo::DiagonalIterator<RANK,V> end  () const;
  //@}
  
protected:
  LazyDensityOperator(const Dimensions& dims) : Base(dims) {}

private:
  virtual const dcomp index(const Idx& i, const Idx& j) const = 0;

};


//@{
  /// Converts the index-tiny of size `RANK` to the arity-dependent indexing type of LazyDensityOperator
template<int RANK>
inline const typename LazyDensityOperator<RANK>::Idx dispatchLDO_index(const TTD_IDXTINY(RANK)& idx) {return idx   ;}

inline const          LazyDensityOperator<1   >::Idx dispatchLDO_index(const TTD_IDXTINY(1   )& idx) {return idx[0];}
//@}


/// Turns the data of a LazyDensityOperator into a real 1D array
/**
 * The result is an array starting with the diagonals (real) of the matrix, optionally followed by the off-diagonals (real & imaginary) of the upper triangle.
 * For this to make sense, the matrix is assumed to be Hermitian.
 * 
 * \param offDiagonals governs whether the upper triangle of the (multi-)matrix is included in the result
 */
template<int RANK>
const DArray<1> deflate(const LazyDensityOperator<RANK>&, bool offDiagonals);


/** \class LazyDensityOperator
 * 
 * \Semantics
 * 
 * - *Unary system*: Assume a mode represented in Fock basis with ladder-operator \f$a\f$. To calculate the quantum expectation value
 * \f[\avr{a^2}=\Tr{a^2\rho}=\sum_i\sqrt{i(i-1)}\,\rho_{i;i-2},\f] one can write the following function:
 *   ~~~
 *   include "LazyDensityOperator.h"
 *   
 *   const dcomp calculateASqr(const LazyDensityOperator<1>& matrix)
 *   {
 *     dcomp res;
 *     for (int i=2; i<matrix.getTotalDimension(); ++i) res+=sqrt(i*(i-1))*matrix(i,i-2);
 *     return res;
 *   }
 *   ~~~
 * 
 * - *Binary system*: Assume two modes represented in Fock bases with ladder-operators \f$a\f$ and \f$b\f$, respectively.
 *   To calculate the quantum expectation value\f[\avr{a^\dag b}=\Tr{a^\dag b\rho}=\sum_{i,j}\sqrt{(i+1)j}\,\rho_{(i,j);(i+1,j-1)},\f]
 *   one can write the following function:
 *   ~~~ 
 *   include "LazyDensityOperator.h"
 * 
 *   const dcomp calculateADaggerB(const LazyDensityOperator<2>& matrix)
 *   {
 *     typedef LazyDensityOperator<2>::Idx Idx;
 *     const LazyDensityOperator<2>::Dimensions dim(matrix.getDimensions());
 * 
 *     dcomp res;
 *     for (int i=0; i<dim[0]-1; ++i) for (int j=1; j<dim[1]; ++j) res+=sqrt((i+1)*j)*matrix(Idx(i,j),Idx(i+1,j-1));
 *     return res;
 *   }
 *   ~~~
 * 
 */

/** \fn template<typename V, typename T, int RANK, typename F> const T partialTrace(const LazyDensityOperator<RANK>&, F function)
 * 
 * \Semantics
 * 
 * The function iterates through all the combinations of the dummy indices, for each slice it takes the value returned by the functor `function`, and accumulates these values.
 * In the following we give some instructive examples of usage.
 * 
 * - *Calculating the full partial density operator of a unary subsystem*: 
 *   The partial density operator of any unary subsystem of a system of any rank may be calculated as follows (e.g. from a StateVector):
 *   ~~~
 *   template<int SUBSYSTEM, int RANK> // for the subsystem indexed by the index SUBSYSTEM
 *   const DensityOperator<1>
 *   partialTraceOfUnarySubsystem(const StateVector<RANK>& psi)
 *   {
 *     return partialTrace<tmptools::Vector<SUBSYSTEM>,DensityOperator<1> >(psi,densityOperatorize<1>);
 *   }
 *   ~~~
 * 
 * - *Calculating correlations for two harmonic-oscillator modes*: Assume having the function `calculateADaggerB` as defined @ LazyDensityOperator.
 *   Then, if these modes are embedded at index positions 3 and 1 (note that the order might be important) in a larger system of arity larger than 3,
 *   we can write the following function:
 *   ~~~
 *   template<int RANK> // RANK>3
 *   const dcomp
 *   calculateADaggerB_atPositions3and1(const LazyDensityOperator<RANK>& matrix)
 *   {
 *     return partialTrace<tmptools::Vector<3,1>,dcomp>(matrix,calculateADaggerB);
 *   }
 *   ~~~
 *   This will calculate \f$\avr{a^\dag b}\f$ (now \f$b\f$ and \f$a\f$ being the ladder operators of the modes at position 1 and 3, respectively)
 *   for the partial density operator of the embedded system (but without explicitly calculating this partial density operator).
 *
 * - *Accumulating an arbitrary ensemble of quantum averages*: An arbitrary ensemble of real or complex numbers can be stored in a DArray<1> or CArray<1>
 *   (or even `std::valarray`s), as these types fulfill all the requirements on type `T`. The former was used to define the abstract interface structure::Averaged
 *   one of its implementations being :class:`Mode`. E.g. a function for a harmonic-oscillator mode calculating \f$\avr{a^\dag a}\f$, \f$\avr{\lp a^\dag a\rp^2}\f$,
 *   and the real and imaginary parts of \f$\avr{a}\f$, may be defined as
 *   ~~~
 *   typedef DArray<1> Averages;
 *
 *   const Averages
 *   calculateModeAverages(const LazyDensityOperator<1>& matrix)
 *   {
 *     Averages averages(4); averages=0;
 *
 *     for (int n=1; n<int(matrix.getDimension()); n++) {
 *       double diag=matrix(n);
 *       averages(0)+=  n*diag;
 *       averages(1)+=n*n*diag;
 *
 *       double sqrtn=sqrt(n);
 *       dcomp offdiag(matrix(n,n-1));
 *       averages(2)+=sqrtn*real(offdiag);
 *       averages(3)+=sqrtn*imag(offdiag);
 *     }
 *
 *     return averages;
 *
 *   }
 *   ~~~
 *   Then the following function will calculate these averages for a mode embedded in a larger system at index position `MODE_POSITION`:
 *   ~~~
 *   template<int RANK, int MODE_POSITION> // for the subsystem indexed by the index MODE_POSITION
 *   const Averages
 *   calculateEmbeddedModeAverages(const StateVector<RANK>& psi)
 *   {
 *     return partialTrace(psi,calculateModeAverages,tmptools::Vector<MODE_POSITION>(),Averages()); 
 *   }
 *   ~~~
 */


} // quantumdata

#endif // QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
