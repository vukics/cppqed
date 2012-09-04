// Randomized GSL implementation

#include "Randomized.h"

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>

#include <boost/make_shared.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


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
    : _ran_gen(gsl_rng_alloc(ran_gen_type)) {gsl_rng_set(_ran_gen,seed);}

  ~RandomizedGSL() {delete _ran_gen;}
  
private:
  double doSample() {return gsl_rng_uniform(_ran_gen);}

  ostream& writeState(ostream&) const;
  istream& readState(istream&);
  ostream& writeImplID(ostream&) const;
  istream& readImplID(istream&) const;

  gsl_rng*const _ran_gen;
  // NEEDS_WORK can this be wrapped into shared_ptr? the same in Evolved.cc?
  static const string _implID;

};


const string RandomizedGSL::_implID("RandomizedGSL");



ostream& RandomizedGSL::writeState(ostream &os) const
{
  using namespace boost::archive::iterators;
  typedef base64_from_binary<  transform_width<  const char *, 6, 8 > > base64_t;

  {
    char *state=static_cast<char*>(gsl_rng_state(_ran_gen));
    copy(base64_t(state),base64_t(state+gsl_rng_size(_ran_gen)),ostream_iterator<char>(os));
  }

  return os << '~' ;
}


istream& RandomizedGSL::readState(istream &is)
{
  using namespace boost::archive::iterators;
  typedef transform_width< binary_from_base64< istream_iterator<char> >, 8, 6 > binary_t;

  {
    char *state=static_cast<char*>(gsl_rng_state(_ran_gen));
    binary_t begin = binary_t(istream_iterator<char>(is));
    // don't increment anymore than necessary
    for ( int size = gsl_rng_size(_ran_gen); --size>0; *state++=static_cast<char>(*(begin++)) );
    *state++ = static_cast<char>(*begin);
  }

  {
    char c = is.get();
    if (c!='~') throw RNGStateParsingException(string("Expected ~, found ")+string(1,c));
  }

  return is;

}


ostream& RandomizedGSL::writeImplID(ostream &os) const
{
  return os << _implID << " ";
}


istream& RandomizedGSL::readImplID(istream &is) const
{
  string id;
  char c;
  is >> id;
  if (id!=_implID) throw RNGStateParsingException(string("Wrong implementation ID, expected "+_implID+", found "+id));
  if (!(is >> c)) throw RNGStateParsingException(string("Unexpected EOF"));
  return is.unget();
}


const Randomized::SmartPtr MakerGSL::operator()(unsigned long seed) const
{
  return boost::make_shared<RandomizedGSL>(seed);
}


}
