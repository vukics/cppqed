// -*- C++ -*-
#ifndef _RANDOMIZED_H
#define _RANDOMIZED_H

#include "RandomizedFwd.h"

#include "ComplexExtensions.h"

#include<boost/shared_ptr.hpp>
#include<boost/utility.hpp>


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
  virtual double operator()() const = 0;

  virtual double gaussian(double sigma=1) const = 0;

  dcomp dcompRan() const {double reRan=(*this)(), imRan=(*this)(); return dcomp(reRan,imRan);}
  // Note that the result of return dcomp((*this)(),(*this)()); is undefined!

};


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


} // randomized

#endif // _RANDOMIZED_H
