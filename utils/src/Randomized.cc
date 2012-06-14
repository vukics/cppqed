// Randomized GSL implementation

#include "Randomized.h"

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


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
    : _ran_gen(gsl_rng_alloc(ran_gen_type)), _implID("RandomizedGSL") {gsl_rng_set(_ran_gen,seed);}

  ~RandomizedGSL() {delete _ran_gen;}
  
private:
  double doSample() {return gsl_rng_uniform(_ran_gen);}
  std::ostream& writeState(std::ostream&) const;
  std::istream& readState(std::istream&);
  std::ostream& writeImplID(std::ostream&) const;
  std::istream& readImplID(std::istream&) const;

  gsl_rng*const _ran_gen;
  // NEEDS_WORK can this be wrapped into shared_ptr? the same in Evolved.cc?
  const std::string _implID;
};

std::ostream& RandomizedGSL::writeState(std::ostream &os) const
{
  using namespace boost::archive::iterators;
  typedef base64_from_binary<  transform_width<  const char *, 6, 8 > > base64_t;
  char *state=static_cast<char*>(gsl_rng_state(_ran_gen));
  int size = gsl_rng_size(_ran_gen);
  std::copy(base64_t(state),base64_t(state+size),std::ostream_iterator<char>(os));
  os << '~';
  return os;
}

std::istream& RandomizedGSL::readState(std::istream &is)
{
  using namespace boost::archive::iterators;
  typedef transform_width< binary_from_base64< std::istream_iterator<char> >, 8, 6 > binary_t;
  char *state=static_cast<char*>(gsl_rng_state(_ran_gen));
  int size = gsl_rng_size(_ran_gen);
  binary_t begin = binary_t(std::istream_iterator<char>(is));
  // don't increment anymore than necessary
  while(--size>0) {
    *state++ = static_cast<char>(*begin);
    ++begin;
  }
  *state++ = static_cast<char>(*begin);
  {
    char c = is.get();
    if (c!='~') throw RNGStateParsingException(std::string("Expected ~, found ")+std::string(1,c));
  }
  return is;
}

std::ostream& RandomizedGSL::writeImplID(std::ostream &os) const
{
  os << _implID << " ";
  return os;
}

std::istream& RandomizedGSL::readImplID(std::istream &is) const
{
  std::string id;
  char c;
  is >> id;
  if (id!=_implID) throw RNGStateParsingException(std::string("Wrong implementation ID, expected "+_implID+", found "+id));
  if (!(is >> c)) throw RNGStateParsingException(std::string("Unexpected EOF"));
  is.unget();
  return is;
}

const Randomized::SmartPtr MakerGSL::operator()(unsigned long seed) const
{
  return Randomized::SmartPtr(new RandomizedGSL(seed));
}

}
