// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines blitzplusplus::TinyOfArrays}
#if !BOOST_PP_IS_ITERATING

#ifndef CPPQEDCORE_UTILS_BLITZTINYOFARRAYS_H_INCLUDED
#define CPPQEDCORE_UTILS_BLITZTINYOFARRAYS_H_INCLUDED

#include "BlitzTinyOfArraysFwd.h"

#include "BlitzArray.h"

#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>

#include <boost/operators.hpp>

namespace blitzplusplus {


struct ShallowCopy {}; // For referencing constructors
struct DeepCopy    {}; // For copying     constructors


template<typename T, int RANK, int LENGTH>
bool operator==(const TinyOfArrays<T,RANK,LENGTH>&, const TinyOfArrays<T,RANK,LENGTH>&);


#define BASE_class blitz::TinyVector<blitz::Array<T,RANK>,LENGTH>

/// Adresses the problem discussed in [this thread](http://sourceforge.net/mailarchive/message.php?msg_id=20017905).
template<typename T, int RANK, int LENGTH>
class TinyOfArrays : public BASE_class, private boost::equality_comparable<TinyOfArrays<T,RANK,LENGTH> >
{
public:
  typedef BASE_class Base;

#undef  BASE_class

  typedef typename Base::T_numtype      T_numtype;
  typedef typename Base::T_vector       T_vector;
  typedef typename Base::T_iterator     T_iterator;
  typedef typename Base::iterator       iterator;
  typedef typename Base::const_iterator const_iterator;
  
  using Base::operator=;

  // In all the constructors the base is implicitly default constructed

  TinyOfArrays()  {}
  ~TinyOfArrays() {}

  TinyOfArrays(const TinyOfArrays&);
  // NEED_TO_UNDERSTAND if this is not declared then OK (how, when the composed constructor doesn't have good semantics?), if this is declared private then compile error with return statements, if it is declared public but not DEFINED then OK.
  // For the solution cf. GotW.#1

  inline TinyOfArrays(ShallowCopy, const TinyOfArrays&);
  inline TinyOfArrays(DeepCopy   , const TinyOfArrays&);

  inline TinyOfArrays(ShallowCopy, const T_numtype& initValue);
  inline TinyOfArrays(DeepCopy   , const T_numtype& initValue);


#define BOOST_PP_ITERATION_LIMITS (2,BLITZ_ARRAY_LARGEST_RANK)
#define BOOST_PP_FILENAME_1 "BlitzTinyOfArrays.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS


};



template<typename T, int RANK, int LENGTH> 
const TinyOfArrays<T,RANK,LENGTH>
negate(const TinyOfArrays<T,RANK,LENGTH>&);



} // blitzplusplus


#endif // CPPQEDCORE_UTILS_BLITZTINYOFARRAYS_H_INCLUDED


#else  // BOOST_PP_IS_ITERATING


#define ITER BOOST_PP_ITERATION()

#define print1(z,m,unused) (*this)(m).reference(x##m);
#define print2(z,m,unused) (*this)(m).reference(x##m.copy());


TinyOfArrays(ShallowCopy, BOOST_PP_ENUM_PARAMS(ITER,const T_numtype& x) )
{
  BOOST_PP_REPEAT(ITER,print1,~);
}

TinyOfArrays(DeepCopy   , BOOST_PP_ENUM_PARAMS(ITER,const T_numtype& x) )
{
  BOOST_PP_REPEAT(ITER,print2,~);
}


#undef print2
#undef print1

#undef ITER


#endif // BOOST_PP_IS_ITERATING
