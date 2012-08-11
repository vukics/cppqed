// -*- C++ -*-
#ifndef _RANDOMIZED_H
#define _RANDOMIZED_H

#include "RandomizedFwd.h"

#include "ComplexExtensions.h"
#include "Exception.h"
#include "Range.h"

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#ifdef USE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#include <sstream>
#endif // USE_BOOST_SERIALIZATION


namespace randomized {

class RNGStateParsingException : public cpputils::TaggedException
{
public:
  RNGStateParsingException(const std::string& tag) : cpputils::TaggedException(tag) {}
};
  
///////////////////////
//
// Randomized interface
//
///////////////////////

class Randomized : private boost::noncopyable
{
public:
  // The logical state of the class is the state of the underlying generator, so that everything that (may) change this state, for example sampling, is logically non-const.

  typedef boost::shared_ptr<Randomized> SmartPtr;

  virtual ~Randomized() {}

  double operator()() {return doSample();};

  const dcomp dcompRan();
  
  friend std::ostream& operator<<(std::ostream&, const Randomized&);
  friend std::istream& operator>>(std::istream&, Randomized&);

private:
  virtual double doSample() = 0;
  virtual std::ostream& writeState(std::ostream&) const = 0;
  virtual std::istream& readState(std::istream&) = 0;
  virtual std::ostream& writeImplID(std::ostream&) const = 0;
  virtual std::istream& readImplID(std::istream&) const = 0;

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

inline
std::ostream& operator<<(std::ostream& os, const Randomized &r)
{
  return r.writeState(r.writeImplID(os));
}

inline
std::istream& operator>>(std::istream& is, Randomized &r)
{
  return r.readState(r.readImplID(is));
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
