// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED

#include "ArrayBase.h" // for ByReference
#include "QuantumDataFwd.h"

#include "DimensionsBookkeeper.h"

#include "BlitzArray.h"
#include "BlitzTinyExtensions.h"
#include "ComplexExtensions.h"
#include "MathExtensions.h"
#include "SliceIterator.tcc"

#include <memory>
#include <numeric>


/// Comprises classes representing the state of composite quantum systems and providing various interfaces to manipulate this data
/** Some of its most important classes fit into a single class rooted in the virtual interface LazyDensityOperator. */
namespace quantumdata {


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
  : public DimensionsBookkeeper<RANK>,
    std::enable_shared_from_this<const LazyDensityOperator<RANK>>
{
public:
  typedef std::shared_ptr<const LazyDensityOperator> Ptr; ///< Many class templates in the framework define shared pointers to their own types, in a template-metafunction like manner
  
  typedef DimensionsBookkeeper<RANK> Base;

  typedef typename Base::Dimensions Dimensions; ///< Inherited from DimensionsBookkeeper

  typedef IdxTiny<RANK> Idx; ///< The type used for indexing the “rows” and the “columns”: a tiny vector of integers (multi-index)

  virtual ~LazyDensityOperator() {}

private:
  class IndexerProxy
  {
  public:
    IndexerProxy(const LazyDensityOperator* ldo, const Idx& firstIndex) : ldo_(ldo), firstIndex_(firstIndex) {}

    const dcomp operator()(const Idx& secondIndex) const {return ldo_->index(firstIndex_,secondIndex);}

    template<typename... SubscriptPack>
    const dcomp operator()(int s0, SubscriptPack... subscriptPack) const
    {
      static_assert( sizeof...(SubscriptPack)==RANK-1 , "Incorrect number of subscripts for LazyDensityOperator::IndexerProxy." );
      return operator()(Idx(s0,subscriptPack...));
    }

    operator double() const {return real(ldo_->index(firstIndex_,firstIndex_));}

  private:
    const LazyDensityOperator*const ldo_;
    const Idx firstIndex_;

  };

public:
  /// Multi-matrix style indexing via Idx type
  const IndexerProxy operator()(const Idx& firstIndex) const {return IndexerProxy(this,firstIndex);}

  /// Multi-matrix style indexing via packs of integers
  /**
   * This allows for the very convenient indexing syntax (e.g. a ternary LazyDensityOperator `matrix` indexed by multi-indeces `i` and `j`):
   *
   *     matrix(i0,i1,i2)(j0,j1,j2)
   *
   * while
   *
   *     matrix(i0,i1,i2)
   *
   * returns a proxy implicitly convertible to a `double`, giving the diagonal element corresponding to the multi-index `i`
   *
   * \note The number of indeces in the multi-index is checked @ compile time.
   *
   */
  template<typename... SubscriptPack>
  const auto operator()(int s0, SubscriptPack... subscriptPack) const
  {
    static_assert( sizeof...(SubscriptPack)==RANK-1 , "Incorrect number of subscripts for LazyDensityOperator." );
    return operator()(Idx(s0,subscriptPack...));
  }

  double trace() const {return trace_v();} ///< Returns the trace (redirected to a pure virtual)
  
protected:
  LazyDensityOperator(const Dimensions& dims) : Base(dims) {}

private:
  virtual const dcomp index(const Idx& i, const Idx& j) const = 0;
  
  virtual double trace_v() const = 0;
  
};


template<typename V, template <int> class MATRIX, int RANK, typename F>
auto partialTrace(const MATRIX<RANK>& matrix, F&& function)
{
  auto begin{cpputils::sliceiterator::begin<V,MATRIX>(matrix)};

  std::cerr<<function(*begin)<<std::endl;
  
  return std::accumulate(++begin,
                         cpputils::sliceiterator::end<V,MATRIX>(matrix),
                         function(*begin),
                         [f{std::move(function)}](const auto& res, const auto& slice){std::cerr<<res<<f(slice)<<std::endl; return res+f(slice);}
                        );
}

struct LazyDensityOperatorPartialTraceNotImplementedException : cpputils::Exception {};

/// The primary tool for performing \ref slicinganldo "slice iterations"
/**
 * On higher levels of the framework (cf. eg. BinarySystem, Composite), this function is used exclusively for performing LazyDensityOperator slice iteration.
 * 
 * \tparamV
 * \tparam F a callable type with signature `T(const LazyDensityOperator<mpl::size<V>::value>&)`, where T is an arithmetic type
 * 
 */
template<typename V, int RANK, typename F>
auto
partialTrace(const LazyDensityOperator<RANK>& matrix, F&& function)
{
  typedef StateVector    <RANK> SV;
  typedef DensityOperator<RANK> DO;
  
  if      (auto psi=dynamic_cast<const SV*>(&matrix) ) return partialTrace<V,StateVector    >(*psi,std::move(function));
  else if (auto rho=dynamic_cast<const DO*>(&matrix) ) return partialTrace<V,DensityOperator>(*rho,std::move(function));
  else throw LazyDensityOperatorPartialTraceNotImplementedException();
}


/// Turns the data of a LazyDensityOperator into a real 1D array
/**
 * The result is an array starting with the diagonals (real) of the matrix, optionally followed by the off-diagonals (real & imaginary) of the upper triangle.
 * For this to make sense, the matrix is assumed to be Hermitian.
 * 
 * \param offDiagonals governs whether the upper triangle of the (multi-)matrix is included in the result
 */
template<int RANK>
const DArray<1> deflate(const LazyDensityOperator<RANK>& matrix, bool offDiagonals)
{
  using mathutils::sqr;
  
  const size_t dim=matrix.getTotalDimension();
  
  DArray<1> res(offDiagonals ? sqr(dim) : dim);
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(matrix.getDimensions()-1,cpputils::mii::begin);
  
  size_t idx=0;

  for (Iterator i(etalon); idx<dim; ++i)
    res(idx++)=matrix(*i);
  
  if (offDiagonals)
    for (Iterator i(etalon); idx<sqr(dim); ++i)
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j) {
        dcomp matrixElement(matrix(*i)(*j));
        res(idx++)=real(matrixElement);
        res(idx++)=imag(matrixElement);
      }

  return res;

}


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
 *     const LazyDensityOperator<2>::Dimensions dim(matrix.getDimensions());
 * 
 *     dcomp res;
 *     for (int i=0; i<dim[0]-1; ++i) for (int j=1; j<dim[1]; ++j) res+=sqrt((i+1)*j)*matrix(i,j)(i+1,j-1);
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
 * - *Accumulating an arbitrary set of quantum averages*: An arbitrary set of real or complex numbers can be stored in a DArray<1> or CArray<1>
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
 * It is implemented in terms of a cpputils::SliceIterator, which is adapted to work with StateVector or DensityOperator (or their non-orthogonal counterparts). 
 * A runtime implementation selection occurs in partialTrace.
 * 
 */

} // quantumdata


#endif // CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
