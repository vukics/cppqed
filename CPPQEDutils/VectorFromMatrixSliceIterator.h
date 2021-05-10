// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Tools for adapting blitzplusplus::basi::Iterator%s to iteration over rows or columns of (multi)matrices}
#ifndef CPPQEDCORE_UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED

#include "BlitzArray.h"
#include "SliceIterator.h"


namespace blitzplusplus {


/// The name of the namespace stands for <strong>V</strong>ector<strong>F</strong>rom<strong>M</strong>atrix<strong>S</strong>lice<strong>I</strong>terator
/**
 * It comprises tools for adapting cppqedutils::SliceIterator%s to iteration over rows or columns of (multi)matrices
 * 
 * \par Semantics
 * Consider the following piece of code:
 * \dontinclude cpcPaperExamples.cc
 * \until <5>&);
 * Assume that `actWithA` is defined in such a way that if `psi` is of type `StateVector<5>`, then `actWithA(psi)` results in \f$A\ket\Psi\f$.
 * \skip composeWithA
 * \until }
 * Then, if `rho` is of type `DensityOperator<5>`, then `composeWithA(rho)` results in \f$A\rho\f$.
 * 
 */
namespace vfmsi {


struct Left  : std::false_type {}; ///< Indicates composition from left
struct Right : std::true_type {}; ///< Indicates composition from right


/// Metafunction returning the appropriate set of retained index positions
/**
 * \tparam RANK the arity of the multivectors resulting from the slicing (the rows/colunms of the multimatrix of total arity `2*RANK`)
 * \tparam S should be either Left or Right
 * 
 * \par Semantics
 * 
 * - `LeftRight<RANK,Left>` is equivalent to a compile-time \link tmptools::Range range\endlink of `<0,1,2,3,...,RANK-1>`
 * - `LeftRight<RANK,Right>` ” `<RANK,RANK+1,...,2*RANK-1>`
 */
template<int RANK, typename S>
using LeftRight = tmptools::Range<tmptools::integral_if_v<S,RANK,0>,
                                  tmptools::integral_if_v<S,2*RANK,RANK>>;


/// Template alias
/**
 * \tparam RANK the arity of the multivectors resulting from the slicing (the rows/colunms of the multimatrix of total arity `2*RANK`)
 * \tparam S should be either Left or Right
 * 
 */
template<int RANK, typename S> using Iterator=cppqedutils::SliceIterator<CArray,RANK,LeftRight<RANK/2,S>>;

/// \name Makers for Iterator
//@{
/// Same as begin but here it returns an Iterator instance
template<typename S, typename A>
auto
begin(const A& array ) {return Iterator<cppqedutils::Rank_v<A>,S>(array,cppqedutils::sliceiterator::Begin{});}

template<typename S, typename A>
auto
end  (const A& array ) {return Iterator<cppqedutils::Rank_v<A>,S>(array,cppqedutils::sliceiterator::End{});}

template<typename S, typename A>
auto
fullRange(const A& array ) {return boost::iterator_range<Iterator<cppqedutils::Rank_v<A>,S> >(begin<S>(array),end<S>(array));}
//@}


} // vfmsi

} // blitzplusplus




#endif // CPPQEDCORE_UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED
