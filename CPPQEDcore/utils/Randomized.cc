// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// Randomized GSL implementation

#include "Randomized.h"

#include <boost/make_shared.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <sstream>
#include <cstdio>


using namespace std;


namespace randomized {


const dcomp Randomized::dcompRan()
{
  // Note that the result of return dcomp((*this)(),(*this)()); is undefined!
  double reRan=doSample(), imRan=doSample();
  return dcomp(reRan,imRan);
}


class RandomizedGSL : public Randomized 
{
public: 
  explicit RandomizedGSL(unsigned long seed, const gsl_rng_type* ran_gen_type=gsl_rng_taus2) 
    : ranGen_(gsl_rng_alloc(ran_gen_type)) {gsl_rng_set(ranGen_.get(),seed);}

  ~RandomizedGSL() {}
  
private:
  typedef boost::shared_ptr<gsl_rng> Impl;
  
  double doSample() {return gsl_rng_uniform(ranGen_.get());}

  const string getState() const
  {
    ostringstream stream;
    stream.write(static_cast<const char*>(gsl_rng_state(ranGen_.get())),gsl_rng_size(ranGen_.get()));
    return stream.str();
  }
  
  void setState(const string& stateIn)
  {
    istringstream stream(stateIn);
    stream.read(static_cast<char*>(gsl_rng_state(ranGen_.get())),gsl_rng_size(ranGen_.get()));
  }

  const string getImplID() const {return "RandomizedGSL";}
  
  const Impl ranGen_;

};


const Randomized::Ptr MakerGSL::operator()(unsigned long seed) const
{
  return boost::make_shared<RandomizedGSL>(seed);
}


} // randomized
