// Randomized GSL implementation

#include "Randomized.h"

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>


namespace randomized {

class RandomizedGSL : public Randomized 
{
public: 
  explicit RandomizedGSL(unsigned long seed, const gsl_rng_type* ran_gen_type=gsl_rng_taus2) 
    : _ran_gen(gsl_rng_alloc(ran_gen_type)) {gsl_rng_set(_ran_gen,seed);}

  ~RandomizedGSL() {delete _ran_gen;}
  
  double operator()() const {return gsl_rng_uniform(_ran_gen);}

  double gaussian(double sigma) const {return gsl_ran_gaussian_ziggurat(_ran_gen,sigma);}

private:
  gsl_rng*const _ran_gen;
  // NEEDS_WORK can this be wrapped into shared_ptr? the same in Evolved.cc?
};


const Randomized::SmartPtr MakerGSL::operator()(unsigned long seed) const
{
  return Randomized::SmartPtr(new RandomizedGSL(seed));
}

}
