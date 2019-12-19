// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPPQEDCORE_QUANTUMDATA_TRANSFORMATION_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_TRANSFORMATION_H_INCLUDED

#include "TransformationFwd.h"

#include "Types.h"

#include "CMatrix.h"
#include "BlitzArraySliceIterator.h"

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/algorithm/transformation/join.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>



namespace quantumdata {


/// Comprises tools for metric transformations of NonOrthogonalStateVector and NonOrthogonalDensityOperator classes.
namespace transformation {


template<int RANK>
struct Identity
{
};




// Special Traits:

// The following is a base class for all Traits that describe such transformations as cannot be factorized.

template<typename TRAFO>
struct ElementaryTraits
{
  typedef typename boost::fusion::vector<TRAFO> TrafoTypes;
  // Simply creates a single-element vector out of TRAFO, so that afterwards it can be treated in a unified way with Composites, e.g. by Compose

  static const TrafoTypes getTrafos(const TRAFO& trafo) {return boost::fusion::make_vector(trafo);}
  // The same at runtime

};



template<int RANK>
struct Traits<Identity<RANK> > : ElementaryTraits<Identity<RANK> >
{
  static const int N_RANK=RANK;

  typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;
  static void transform(Identity<RANK>, const StateVectorLow& in, StateVectorLow& out) {out=in;}

};



// Elementary examples are ...
// ... Transformations described by (multi)matrices:

template<int TWO_TIMES_RANK>
struct Traits<CArray<TWO_TIMES_RANK> > : ElementaryTraits<CArray<TWO_TIMES_RANK> >
{
  static const int N_RANK=TWO_TIMES_RANK/2;

  typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;
  static void transform(const CArray<TWO_TIMES_RANK>& trafo, const StateVectorLow& in, StateVectorLow& out);

};


// ... Transformations described by functions of the appropriate signature:

/** \todo For quantumdata::transformation::Traits default implementation should pertain to a function(object). Probably not very easy to implement, though… */
template<int RANK>
struct Traits< void(*)(const CArray<RANK>&, CArray<RANK>&) > 
  : ElementaryTraits< void(*)(const CArray<RANK>&, CArray<RANK>&) >
{
  static const int N_RANK=RANK;

  typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;

  typedef void(*TRAFO)(const StateVectorLow&, StateVectorLow&);

  static void transform(const TRAFO trafo, const StateVectorLow& in, StateVectorLow& out) {trafo(in,out);}

};


/*
// ... Transformations described by functors of the appropriate signature:

template<int RANK>
struct Traits<boost::function<void(const CArray<RANK>&, CArray<RANK>&)> > 
  : ElementaryTraits<boost::function<void(const CArray<RANK>&, CArray<RANK>&)> >
{
  static const int N_RANK=RANK;

  typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;

  typedef void(*TRAFO)(const StateVectorLow&, StateVectorLow&);

  typedef boost::function<void(const CArray<RANK>&, CArray<RANK>&)> TRAFO;

  static void transform(TRAFO trafo, const StateVectorLow& in, StateVectorLow& out) {trafo(in,out);}

};
*/

} // transformation

} // quantumdata

#endif // CPPQEDCORE_QUANTUMDATA_TRANSFORMATION_H_INCLUDED

