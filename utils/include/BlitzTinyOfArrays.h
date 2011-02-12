// -*- C++ -*-
#ifndef _BLITZ_TINY_OF_ARRAYS_H
#define _BLITZ_TINY_OF_ARRAYS_H

#include "BlitzTinyOfArraysFwd.h"

#include<blitz/array.h>

#include<boost/preprocessor/iteration/iterate.hpp>
#include<boost/preprocessor/repetition.hpp>
#include<boost/preprocessor/arithmetic/mul.hpp>


namespace blitzplusplus {


struct TOA_ShallowCopy {}; // For referencing constructors
struct TOA_DeepCopy    {}; // For copying     constructors


#define BASE_class blitz::TinyVector<blitz::Array<T,RANK>,LENGTH>

template<typename T, int RANK, int LENGTH>
class TinyOfArrays : public BASE_class
{
public:
  typedef BASE_class Base;

#undef  BASE_class

  typedef typename Base::T_numtype       T_numtype;
  typedef typename Base::T_vector        T_vector;
  typedef typename Base::T_iterator      T_iterator;
  typedef typename Base::T_constIterator T_constIterator;
  typedef typename Base::iterator        iterator;
  typedef typename Base::const_iterator  const_iterator;
  
  using Base::numElements; using Base::operator=;


  // In all the constructors the base is implicitly default constructed

  TinyOfArrays()  {}
  ~TinyOfArrays() {}

  TinyOfArrays(const TinyOfArrays&);
  // NEED_TO_UNDERSTAND if this is not declared then OK (how, when the
  // composed constructor doesn't have good semantics?), if this is
  // declared private then compile error with return statements, if it
  // is declared public but not DEFINED then OK. For the solution
  // cf. GotW.#1

  inline TinyOfArrays(TOA_ShallowCopy, const TinyOfArrays&);
  inline TinyOfArrays(TOA_DeepCopy   , const TinyOfArrays&);

  inline TinyOfArrays(TOA_ShallowCopy, const T_numtype& initValue);
  inline TinyOfArrays(TOA_DeepCopy   , const T_numtype& initValue);


#define BOOST_PP_ITERATION_LIMITS (2,11)
#define BOOST_PP_FILENAME_1 "details/BlitzTinyOfArraysSpecialization.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS


};



template<typename T, int RANK, int LENGTH> 
const TinyOfArrays<T,RANK,LENGTH>
negate(const TinyOfArrays<T,RANK,LENGTH>&);



} // blitzplusplus


#include "impl/BlitzTinyOfArrays.tcc"


#endif // _BLITZ_TINY_OF_ARRAYS_H
