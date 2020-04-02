// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Tools for adapting blitzplusplus::basi::Iterator%s to iteration over rows or columns of (multi)matrices}
// -*- C++ -*-
#ifndef CPPQEDCORE_UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED

#include "VectorFromMatrixSliceIteratorFwd.h"

#include "BlitzArraySliceIterator.h"


namespace blitzplusplus {


/// The name of the namespace stands for <strong>V</strong>ector<strong>F</strong>rom<strong>M</strong>atrix<strong>S</strong>lice<strong>I</strong>terator
/**
 * It comprises tools for adapting basi::Iterator%s to iteration over rows or columns of (multi)matrices
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


struct Left  : boost::mpl::false_ {}; ///< Indicates composition from left
struct Right : boost::mpl:: true_ {}; ///< Indicates composition from right


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
struct LeftRight 
  : tmptools::Range<RANK,
                    boost::mpl::if_<S,
                                    boost::mpl::int_<RANK>,
                                    boost::mpl::int_<0>
                                    >::type::value> {};


/// Template alias
/**
 * \tparam RANK the arity of the multivectors resulting from the slicing (the rows/colunms of the multimatrix of total arity `2*RANK`)
 * \tparam S should be either Left or Right
 * \tparam IS_CONST governs constness
 * 
 */
template<int RANK, typename S, bool IS_CONST> using Iterator=basi::Iterator<RANK,LeftRight<RANK/2,S>,IS_CONST>;

/// \name Makers for Iterator
//@{
/// Same as begin but here it returns an Iterator instance
template<typename S, typename A>
const Iterator<Rank<A>::value,S,false>
begin(     A& array );

template<typename S, typename A>
const Iterator<Rank<A>::value,S,false>
end  (     A& array );

template<typename S, typename A>
const Iterator<Rank<A>::value,S,true>
begin(const A& array );

template<typename S, typename A>
const Iterator<Rank<A>::value,S,true>
end  (const A& array );

template<typename S, typename A>
const boost::iterator_range<Iterator<Rank<A>::value,S,true> >
fullRange(const A& array );

template<typename S, typename A>
const boost::iterator_range<Iterator<Rank<A>::value,S,false> >
fullRange(      A& array );
//@}


#define NS_NAME vfmsi
#define RETURN_type1(IS_CONST) Iterator<Rank<A>::value,V_S,IS_CONST>
#define ADDITIONAL_PARAMETER
#define ADDITIONAL_ARGUMENT

#include "details_BlitzArraySliceIteratorReentrant.h"


} // vfmsi


/** \cond SPECIALIZATION */

namespace basi {

template<int RANK, typename S> struct ConsistencyChecker<RANK,blitzplusplus::vfmsi::LeftRight<RANK/2,S> > {};

} // basi

/** \endcond */

} // blitzplusplus


#endif // CPPQEDCORE_UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED
