// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Random & related}
#ifndef CPPQEDCORE_UTILS_RANDOM_H_REENTRANT

#ifndef CPPQEDCORE_UTILS_RANDOM_H_INCLUDED
#define CPPQEDCORE_UTILS_RANDOM_H_INCLUDED

#ifdef CPPQED_HAS_GSL
#include <gsl/gsl_rng.h>
#endif // CPPQED_HAS_GSL

#include <boost/range/algorithm/generate.hpp>

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>
#endif // BZ_HAVE_BOOST_SERIALIZATION

#include <random>
#include <memory>
#include <stdexcept>
#include <sstream>


namespace cpputils {

#ifdef CPPQED_HAS_GSL

/// Wraps GSL random number generators into the concept [UniformRandomBitGenerator](http://en.cppreference.com/w/cpp/named_req/UniformRandomBitGenerator)
class GSL_RandomEngine 
{
public:
  using result_type=unsigned long;
  
  explicit GSL_RandomEngine(result_type s/*, const gsl_rng_type* ran_gen_type=gsl_rng_taus2*/);

  GSL_RandomEngine(const GSL_RandomEngine&);

  void seed(result_type value);

  result_type operator()();
  
  result_type min() const;
  result_type max() const;

  void write(std::ostream& os) const;

  void read(std::istream& is);
  
private:
  const std::shared_ptr<gsl_rng> ranGen_;
  const std::string name_;
  
};


inline std::ostream& operator<<(std::ostream& os, const GSL_RandomEngine& r) {r.write(os); return os;}
inline std::istream& operator>>(std::istream& is, GSL_RandomEngine& r) {r.read(is); return is;}


#define CPPQEDCORE_UTILS_RANDOM_H_REENTRANT
#define CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE GSL_RandomEngine
#define CPPQEDCORE_UTILS_RANDOM_H_STRING "GSL_RandomEngine state"

#include "Random.h"

} // cpputils

#endif // CPPQED_HAS_GSL


namespace cpputils {

template<typename Distribution, typename RandomEngine, typename Array, typename... DCP>
RandomEngine& fillWithRandom(Array a, RandomEngine& re, DCP&&... dcp)
{
  Distribution d{std::forward<DCP>(dcp)...};
  std::generate(a.begin(),a.end(),[&]() {return d(re);});
  return re;
}


template<typename Distribution, typename RandomEngine, typename Array, typename... DCP>
void fillWithRandom(Array a, unsigned long seed, DCP&&... dcp)
{
  RandomEngine re{seed};
  fillWithRandom<Distribution>(a,re,std::forward<DCP>(dcp)...);
}

} // cpputils


namespace boost { namespace serialization {

#define CPPQEDCORE_UTILS_RANDOM_H_REENTRANT
#define CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE std::mt19937_64
#define CPPQEDCORE_UTILS_RANDOM_H_STRING "mersenne_twister_engine state"

#include "Random.h"

} } // boost::serialization

#endif // CPPQEDCORE_UTILS_RANDOM_H_INCLUDED


#else // CPPQEDCORE_UTILS_RANDOM_H_REENTRANT


template<typename Ar>
void load(Ar& ar, CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE& mt, unsigned) {
  std::string text;
  ar & text;
  std::istringstream iss(text);

  if (!(iss >> mt))
    throw std::invalid_argument(CPPQEDCORE_UTILS_RANDOM_H_STRING);
}

template<typename Ar>
void save(Ar& ar, CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE const& mt, unsigned) {
  std::ostringstream oss;
  if (!(oss << mt))
    throw std::invalid_argument(CPPQEDCORE_UTILS_RANDOM_H_STRING);
  std::string text = oss.str();
  ar & text;
}

template<typename Ar>
void serialize(Ar& ar, CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE& mt, unsigned version) {
  if (typename Ar::is_saving())
    save(ar, mt, version);
  else
    load(ar, mt, version);
}

#undef CPPQEDCORE_UTILS_RANDOM_H_STRING
#undef CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE
#undef CPPQEDCORE_UTILS_RANDOM_H_REENTRANT

#endif // CPPQEDCORE_UTILS_RANDOM_H_REENTRANT
