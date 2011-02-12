// -*- C++ -*-
#ifndef _DUAL_TRANSFORMATION_H
#define _DUAL_TRANSFORMATION_H

#include "TransformationFwd.h"

#include "Types.h"

#include "Algorithm.h"
#include "CMatrix.h"
#include "BlitzArraySliceIterator.h"
#include "Range.h"

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/algorithm/transformation/join.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>



namespace quantumdata {

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
struct Traits<TTD_CARRAY(TWO_TIMES_RANK)> : ElementaryTraits<TTD_CARRAY(TWO_TIMES_RANK)>
{
  static const int N_RANK=TWO_TIMES_RANK/2;

  typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;
  static void transform(const TTD_CARRAY(TWO_TIMES_RANK)& trafo, const StateVectorLow& in, StateVectorLow& out);

};


// ... Transformations described by functions of the appropriate signature:

template<int RANK>
struct Traits< void(*)(const TTD_CARRAY(RANK)&, TTD_CARRAY(RANK)&) > 
  : ElementaryTraits< void(*)(const TTD_CARRAY(RANK)&, TTD_CARRAY(RANK)&) >
{
  static const int N_RANK=RANK;

  typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;

  typedef void(*TRAFO)(const StateVectorLow&, StateVectorLow&);

  static void transform(const TRAFO trafo, const StateVectorLow& in, StateVectorLow& out) {trafo(in,out);}

};


/*
// ... Transformations described by functors of the appropriate signature:

template<int RANK>
struct Traits<boost::function<void(const TTD_CARRAY(RANK)&, TTD_CARRAY(RANK)&)> > 
  : ElementaryTraits<boost::function<void(const TTD_CARRAY(RANK)&, TTD_CARRAY(RANK)&)> >
{
  static const int N_RANK=RANK;

  typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;

  typedef void(*TRAFO)(const StateVectorLow&, StateVectorLow&);

  typedef boost::function<void(const TTD_CARRAY(RANK)&, TTD_CARRAY(RANK)&)> TRAFO;

  static void transform(TRAFO trafo, const StateVectorLow& in, StateVectorLow& out) {trafo(in,out);}

};
*/

} // transformation

} // quantumdata

#include "impl/Transformation.tcc"

#endif // _DUAL_TRANSFORMATION_H

