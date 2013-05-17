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
/**
 * Some of its most important classes fit into a single class rooted in the virtual interface LazyDensityOperator.
 */
namespace quantumdata {



template<typename V, typename T, int RANK, typename F>
const T
partialTrace(const LazyDensityOperator<RANK>&, F);
// F must be an object callable with the signature:
// const T(const typename ldo::DiagonalIterator<RANK,V>::value_type&)


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
  /**
   * Just an integer (index) if `RANK=1`, otherwise a tiny vector of integers (multi-index).
   */
  typedef typename mpl::if_c<(RANK==1),int,TTD_IDXTINY(RANK)>::type Idx; 

  virtual ~LazyDensityOperator() {}

  const dcomp operator()(const Idx& i, const Idx& j) const {return index(i,j);} ///< The indexing function calling a purely virtual function according to the non-virtual interface idiom

  double operator()(const Idx& i) const {return real((*this)(i,i));} ///< An inline function for conveniently addressing the diagonal elements

 
  //@{
    /// Return the DiagonalIterator corresponding to the beginning/end of the sequence of slices defined by `V` 
    /**
     * (cf. the section on slicing a LazyDensityOperator).
     * \tparam V Compile-time vector holding the *retained index positions*.
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


} // quantumdata

#endif // QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
