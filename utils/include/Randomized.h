// -*- C++ -*-
#ifndef UTILS_INCLUDE_RANDOMIZED_H_INCLUDED
#define UTILS_INCLUDE_RANDOMIZED_H_INCLUDED

#include "RandomizedFwd.h"

#include "ComplexExtensions.h"
#include "Exception.h"
#include "Range.h"

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>
#endif // DO_NOT_USE_BOOST_SERIALIZATION


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

  typedef boost::shared_ptr<Randomized> Ptr;

  virtual ~Randomized() {}

  double operator()() {return doSample();};

  const dcomp dcompRan();
  
private:
  virtual double doSample() = 0;

#ifndef DO_NOT_USE_BOOST_SERIALIZATION

  friend class boost::serialization::access;

  template<class Archive>
  void save(Archive& ar, const unsigned int /* version */) const
  {
    const std::string state(getState()), id(getImplID());
    ar  & state & id;
  }
  
  template<class Archive>
  void load(Archive& ar, const unsigned int /* version */)
  {
    std::string state, id;
    ar & state & id;
    if (id!=getImplID()) throw RNGStateParsingException(std::string("Wrong implementation ID, expected "+id+", found "+getImplID()));
    setState(state);
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()
  
#endif // DO_NOT_USE_BOOST_SERIALIZATION
  
  virtual const std::string getState() const = 0;
  virtual void setState(const std::string&) = 0;
  
  virtual const std::string getImplID() const = 0;

};


template<typename D>
inline
const D sample(Randomized::Ptr ran);

template<>
inline
const double sample<double>(Randomized::Ptr ran)
{
  return ran->operator()();
}

template<>
inline
const dcomp  sample<dcomp >(Randomized::Ptr ran)
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
  virtual const Randomized::Ptr operator()(unsigned long seed) const = 0;

  virtual ~Maker() {}
  
};


class MakerGSL : public Maker
{
public:
  const Randomized::Ptr operator()(unsigned long seed) const;
  
};



////////////////////////////
//
// generating a random array
//
////////////////////////////

template<typename A>
const Randomized::Ptr fillWithRandom(A& data, Randomized::Ptr ran)
{
  boost::generate(data,boost::bind(sample<typename A::T_numtype>,ran));
  return ran;
}


template<typename A>
const Randomized::Ptr fillWithRandom(A& data, unsigned long seed=1001ul, const Maker& maker=MakerGSL())
{
  Randomized::Ptr ran(maker(seed));
  return fillWithRandom(data,ran);
}


} // randomized

#endif // UTILS_INCLUDE_RANDOMIZED_H_INCLUDED
