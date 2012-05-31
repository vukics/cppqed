// Randomized GSL implementation

#include "Randomized.h"

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>


namespace randomized {


const dcomp Randomized::dcompRan() const
{
  // Note that the result of return dcomp((*this)(),(*this)()); is undefined!
  double reRan=doSample(), imRan=doSample();
  return dcomp(reRan,imRan);
}


class RandomizedGSL : public Randomized 
{
public: 
  explicit RandomizedGSL(unsigned long seed, const gsl_rng_type* ran_gen_type=gsl_rng_taus2) 
    : _ran_gen(gsl_rng_alloc(ran_gen_type)) {gsl_rng_set(_ran_gen,seed);}

  ~RandomizedGSL() {delete _ran_gen;}
  
private:
  double doSample() const {return gsl_rng_uniform(_ran_gen);}
  std::ostream& writeState(std::ostream&) const;
  std::istream& readState(std::istream&) const;

  gsl_rng*const _ran_gen;
  // NEEDS_WORK can this be wrapped into shared_ptr? the same in Evolved.cc?
};

std::ostream& RandomizedGSL::writeState(std::ostream& os) const
{
  void *state=gsl_rng_state(_ran_gen);
  os.write((char*) state,gsl_rng_size(_ran_gen));
  return os;
}

std::istream& RandomizedGSL::readState(std::istream& is) const
{
  void *state=gsl_rng_state(_ran_gen);
  is.read((char*) state,gsl_rng_size(_ran_gen));
  return is;
}



const Randomized::SmartPtr MakerGSL::operator()(unsigned long seed) const
{
  return Randomized::SmartPtr(new RandomizedGSL(seed));
}

}
