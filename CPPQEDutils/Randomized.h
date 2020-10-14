// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Randomized & related}
#ifndef CPPQEDCORE_UTILS_RANDOMIZED_H_INCLUDED
#define CPPQEDCORE_UTILS_RANDOMIZED_H_INCLUDED

#include "ArrayTraits.h"
#include "ComplexExtensions.h"
#include "Exception.h"

#include <boost/range/algorithm/generate.hpp>

#include <boost/bind.hpp>
#include <boost/utility.hpp>

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>
#endif // BZ_HAVE_BOOST_SERIALIZATION


/// the randomized-bundle
namespace randomized {

class RNGStateParsingException : public cpputils::TaggedException
{
public:
  RNGStateParsingException(const std::string& tag) : cpputils::TaggedException(tag) {}
};


/// A common interface for random-number generators
/**
 * The class can serialize the state of the generator allowing for restoration from an archive (cf. \refBoost{Boost.Serialization,serialization})
 * 
 * \note The logical state of the class is the state of the underlying generator, so that everything that (may) change this state, for example sampling, is logically non-const.
 */
class Randomized : private boost::noncopyable
{
public:
  typedef std::shared_ptr<Randomized> Ptr;

  virtual ~Randomized() {}

  double operator()() {return doSample();} ///< sampling of uniform distribution over the interval [0:1)

  const dcomp dcompRan(); ///< sampling of a uniform distribution over unit square on the complex plane
  
private:
  virtual double doSample() = 0;

#ifdef BZ_HAVE_BOOST_SERIALIZATION

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
    if (id!=getImplID()) throw RNGStateParsingException("Wrong implementation ID, expected "+id+", found "+getImplID());
    setState(state);
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()
  
#endif // BZ_HAVE_BOOST_SERIALIZATION
  
  virtual const std::string getState() const = 0;
  virtual void setState(const std::string&) = 0;
  
  virtual const std::string getImplID() const = 0;

};


/// \related Randomized
template<typename D>
inline
const D sample(Randomized::Ptr ran);

/// \related Randomized
template<>
inline
const double sample<double>(Randomized::Ptr ran)
{
  return ran->operator()();
}

/// \related Randomized
template<>
inline
const dcomp  sample<dcomp >(Randomized::Ptr ran)
{
  return ran->dcompRan();
}


/// Factory class for Randomized types
class Maker
{
public:
  virtual const Randomized::Ptr operator()(unsigned long seed) const = 0;

  virtual ~Maker() {}
  
};


/// Implements Maker by returning a class implementing the Randomized interface by [GSL](http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Distributions.html#Random-Number-Distributions)
class MakerGSL : public Maker
{
public:
  const Randomized::Ptr operator()(unsigned long seed) const;
  
};



/// Fills an array with random data taking a Randomized as parameter
template<typename A>
const Randomized::Ptr fillWithRandom(A& data, Randomized::Ptr ran)
{
  boost::generate(data,boost::bind(sample<typename cpputils::ElementType<A>::type>,ran));
  return ran;
}


/// Fills an array with random data creating a Randomized internally
template<typename A>
const Randomized::Ptr fillWithRandom(A& data, unsigned long seed=1001ul, const Maker& maker=MakerGSL())
{
  Randomized::Ptr ran(maker(seed));
  return fillWithRandom(data,ran);
}


} // randomized

#endif // CPPQEDCORE_UTILS_RANDOMIZED_H_INCLUDED
