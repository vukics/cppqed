// -*- C++ -*-
#ifndef _RANDOMIZED_H
#define _RANDOMIZED_H

#include "RandomizedFwd.h"

#include "ComplexExtensions.h"
#include "Range.h"

#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>


namespace randomized {

///////////////////////
//
// Randomized interface
//
///////////////////////

class Randomized : private boost::noncopyable
{
public:
  typedef boost::shared_ptr<const Randomized> SmartPtr;

  virtual ~Randomized() {}

  double operator()() const {return doSample();};

  const dcomp dcompRan() const;

private:
  virtual double doSample() const = 0;

};


template<typename D>
inline
const D sample(Randomized::SmartPtr ran);

template<>
inline
const double sample<double>(Randomized::SmartPtr ran)
{
  return ran->operator()();
}

template<>
inline
const dcomp  sample<dcomp >(Randomized::SmartPtr ran)
{
  return ran->dcompRan();
}



////////////////
//
// factory class
//
////////////////

class Maker
{
public:
  virtual const Randomized::SmartPtr operator()(unsigned long seed) const = 0;

  virtual ~Maker() {}
  
};


class MakerGSL : public Maker
{
public:
  const Randomized::SmartPtr operator()(unsigned long seed) const;
  
};



////////////////////////////
//
// generating a random array
//
////////////////////////////

template<typename A>
const Randomized::SmartPtr fillWithRandom(A& data, Randomized::SmartPtr ran)
{
  boost::generate(data,boost::bind(sample<typename A::T_numtype>,ran));
  return ran;
}


template<typename A>
const Randomized::SmartPtr fillWithRandom(A& data, unsigned long seed=1001ul, const Maker& maker=MakerGSL())
{
  Randomized::SmartPtr ran(maker(seed));
  return fillWithRandom(data,ran);
}


} // randomized

#endif // _RANDOMIZED_H
