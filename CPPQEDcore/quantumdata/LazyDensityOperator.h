// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "DensityOperator.h"


namespace quantumdata {


/// Common interface for calculating quantum averages from StateVector, DensityOperator (and potentially their non-orthogonal counterparts)
/**
 * Relies on the formula
 * \f[\avr{A}=\Tr{A\rho}\f]
 * (\f$A\f$ is an observable, and \f$\rho\f$ the density operator of the system).
 * 
 * \f$\rho\f$ is simply the dyad of the state vector in the pure-state case, which is however treated lazily (usually not all the matrix elements of \f$\rho\f$ are needed for calculating the average, 
 * plus for large dimensionality, we may not be able to store the full \f$\ket\Psi\bra\Psi\f$).
 * 
 * The solution adopted for this problem in the framework is represented by LazyDensityOperator, which provides a common interface 
 * for all the four cases , to calculate quantum averages from their data.
 * The “laziness” amounts to that in the case of state vectors, only those elements of the density operator are calculated that are actually asked for.
 * 
 * Immutable class. Should be passed by value.
 */


/// The lazy_density_operator-style indexing should be defined in the language of MultiArrayConstView
/**
 * If we do it as DensityOperatorConstView, then template-argument-deduction doesn’t work very well
 * TODO: (re)evaluate the alternative design of making StateVectorConstView & DensityOperatorConstView optional bases of MultiArray
 *   In that case, lazy_density_operator-style indexing could be defined as member functions
 */
inline dcomp _(MultiArrayConstView<dcomp,1> psi, size_t i, size_t j) {return psi(i)*std::conj(psi(j));}

inline double _(MultiArrayConstView<dcomp,1> psi, size_t i) {return sqrAbs(psi(i));}

inline dcomp _(MultiArrayConstView<dcomp,2> rho, size_t i, size_t j) {return rho(i,j);}

inline double _(MultiArrayConstView<dcomp,2> rho, size_t i) {return real(rho(i,i));}

template <size_t RANK>
dcomp _(MultiArrayConstView<dcomp,RANK> psi, Extents<RANK> i, Extents<RANK> j)
{
  return psi(i)*std::conj(psi(j));
}

template <size_t RANK>
double _(MultiArrayConstView<dcomp,RANK> psi, Extents<RANK> i)
{
  return sqrAbs(psi(i));
}

template <size_t TWO_TIMES_RANK>
dcomp _(MultiArrayConstView<dcomp,TWO_TIMES_RANK> rho,
        Extents<TWO_TIMES_RANK/2> i,
        Extents<TWO_TIMES_RANK/2> j) requires (TWO_TIMES_RANK%2==0)
{
  return rho(concatenate(i,j));
}

template <size_t TWO_TIMES_RANK>
double _(MultiArrayConstView<dcomp,TWO_TIMES_RANK> rho,
         Extents<TWO_TIMES_RANK/2> i) requires (TWO_TIMES_RANK%2==0)
{
  return real(rho(concatenate(i,i)));
}



template <typename T, size_t RANK>
concept lazy_density_operator = requires (T&& t, Extents<RANK> i, Extents<RANK> j)
{
  {_(t,i,j)} -> std::convertible_to<dcomp>;
  {_(t,i)} -> std::convertible_to<double>;
};



/// The primary tool for performing slice iteration of LazyDensityOperator#s
/**
 * On higher levels of the framework (BinarySystem, Composite), this function is used exclusively for performing LazyDensityOperator slice iteration.
 *
 * `function` is a callable with signature `T(LazyDensityOperator<hana::size(retainedAxes)>)`, where T is an arithmetic type
 */
template<auto retainedAxes, size_t RANK>
auto partialTrace(lazy_density_operator<RANK> auto matrix, const std::vector<size_t>& offsets, auto&& function, auto&& plus)
{
  // static constexpr auto extendedAxes{hana::fold(retainedAxes,retainedAxes, [] (const auto& state, auto element) {return hana::append(state,element+RANK);})};
  using LDO = decltype(matrix);
 
  const auto iterateSlices{ [&] (const auto& sr) {
    return std::ranges::fold_left_first(sr | std::views::transform(function), plus ).value();
  }};

  if constexpr (std::same_as<LDO,StateVectorConstView<RANK>>) {
    auto sr{sliceRange<retainedAxes>(matrix,offsets)};
    return iterateSlices(sr);
  }
  else if constexpr (std::same_as<LDO,DensityOperatorConstView<RANK>>) {
    static constexpr auto extendedAxes{hana::concat(retainedAxes,
                                                    hana::transform(retainedAxes, [] (const auto& e) {return e+RANK;} ))};
    auto diagonalOffsets{ offsets | std::views::transform( [&] (size_t v) {return v * ( std::lround(std::sqrt(matrix.dataView.size())) + 1 ) ; } ) };
    auto sr{sliceRange<extendedAxes>(matrix,diagonalOffsets)};
    return iterateSlices(sr);
  }
  else static_assert(always_false_v<LDO>, "unknown type input in partialTrace");
 
}


template<auto retainedAxes, size_t RANK>
auto partialTrace(lazy_density_operator<RANK> auto matrix, const std::vector<size_t>& offsets, auto&& function)
{
  return partialTrace<retainedAxes,RANK>(matrix,offsets,std::forward<decltype(function)>(function),std::plus{});
}


/// Turns the data of a LazyDensityOperator into a real 1D array
/**
 * The result is an array starting with the diagonals (real) of the matrix, optionally followed by the off-diagonals (real & imaginary) of the upper triangle.
 * For this to make sense, the matrix is assumed to be Hermitian.
 * 
 * \param offDiagonals governs whether the upper triangle of the (multi-)matrix is included in the result
 *//*
template<size_t RANK>
const DArray<1> deflate(const LazyDensityOperator<RANK>& matrix, bool offDiagonals)
{
  using cppqedutils::sqr;
  
  const size_t dim=matrix.getTotalDimension();
  
  DArray<1> res(offDiagonals ? sqr(dim) : dim);
  
  typedef cppqedutils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(matrix.getDimensions()-1,cppqedutils::mii::begin);
  
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
*/

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

/** \fn template<typename V, typename T, size_t RANK, typename F> const T partialTrace(const LazyDensityOperator<RANK>&, F function)
 * 
 * \Semantics
 * 
 * The function iterates through all the combinations of the dummy indices, for each slice it takes the value returned by the functor `function`, and accumulates these values.
 * In the following we give some instructive examples of usage.
 * 
 * - *Calculating the full partial density operator of a unary subsystem*: 
 *   The partial density operator of any unary subsystem of a system of any rank may be calculated as follows (e.g. from a StateVector):
 *   ~~~
 *   template<int SUBSYSTEM, size_t RANK> // for the subsystem indexed by the index SUBSYSTEM
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
 *   template<size_t RANK> // RANK>3
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
 *   template<size_t RANK, int MODE_POSITION> // for the subsystem indexed by the index MODE_POSITION
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
 * It is implemented in terms of a cppqedutils::SliceIterator, which is adapted to work with StateVector or DensityOperator (or their non-orthogonal counterparts). 
 * A runtime implementation selection occurs in partialTrace.
 * 
 */

} // quantumdata


