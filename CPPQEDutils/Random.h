// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Random & related}
#ifndef CPPQEDCORE_UTILS_RANDOM_H_REENTRANT

#ifndef CPPQEDCORE_UTILS_RANDOM_H_INCLUDED
#define CPPQEDCORE_UTILS_RANDOM_H_INCLUDED

#include "Pars.h"

#include "pcg_random.hpp"
#include "XoshiroCpp.hpp"

#ifdef CPPQED_HAS_GSL
#include <gsl/gsl_rng.h>
#endif // CPPQED_HAS_GSL

#include <boost/range/algorithm/generate.hpp>

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/string.hpp>
#include <boost/serialization/array.hpp>
#endif // BZ_HAVE_BOOST_SERIALIZATION

#include <memory>
#include <random>
#include <stdexcept>
#include <sstream>


namespace randomutils {

template<typename Engine, typename BASE>
struct Pars;


#ifdef CPPQED_HAS_GSL

/// Wraps GSL random number generators into the concept [UniformRandomBitGenerator](http://en.cppreference.com/w/cpp/named_req/UniformRandomBitGenerator)
class GSL_Engine 
{
public:
  using result_type=unsigned long;
  
  explicit GSL_Engine(result_type s/*, const gsl_rng_type* ran_gen_type=gsl_rng_taus2*/);

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


inline std::ostream& operator<<(std::ostream& os, const GSL_Engine& r) {r.write(os); return os;}
inline std::istream& operator>>(std::istream& is, GSL_Engine& r) {r.read(is); return is;}


#define CPPQEDCORE_UTILS_RANDOM_H_REENTRANT
#define CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE GSL_Engine
#define CPPQEDCORE_UTILS_RANDOM_H_STRING "GSL_Engine state"

#include "Random.h"

} // randomutils

#endif // CPPQED_HAS_GSL


namespace randomutils {

template<typename Distribution, typename Engine, typename Array, typename... DCP>
Engine& fill(Array a, Engine& re, DCP&&... dcp)
{
  Distribution d{std::forward<DCP>(dcp)...};
  std::generate(a.begin(),a.end(),[&]() {return d(re);});
  return re;
}


template<typename Distribution, typename Engine, typename Array, typename... DCP>
void fill(Array a, unsigned long seed, DCP&&... dcp)
{
  Engine re{seed};
  fill<Distribution>(a,re,std::forward<DCP>(dcp)...);
}


template<typename Engine>
constexpr auto EngineID_v = std::nullopt;

template<>
inline const std::string EngineID_v<GSL_Engine> = "GSL_Taus2";

template<>
inline const std::string EngineID_v<std::mt19937_64> = "STD_MT19937_64";

template<>
inline const std::string EngineID_v<pcg64> = "PCG64";

template<>
inline const std::string EngineID_v<XoshiroCpp::Xoshiro256PlusPlus> = "Xoshiro256pp";


template<typename Engine>
struct EngineWithParameters
{
  EngineWithParameters(unsigned long s, unsigned long p) : engine{s,p}, seed{s}, prngStream{p} {}
  
  std::ostream& stream(std::ostream& os) const {return os << "Random engine: "<<EngineID_v<Engine><<". Parameters: seed="<<seed<<", streamOrdo=" << prngStream;}

  Engine engine;
  const unsigned long seed, prngStream;

};


// template<>
// inline std::list<pcg64> independentPRNG_Streams<pcg64>(size_t n, pcg64::result_type seed, pcg64::result_type ordoStream)
// {
//   std::list<pcg64> res;
//   for (size_t i=0; i<n; ++i) res.emplace_back(seed,1001+ordoStream+i);
//   return res;
// }

template<typename BASE>
struct Pars<pcg64,BASE> : BASE
{
  unsigned long &seed, &prngStream; ///< PCG64 generator seed and ordinal of independent stream

  Pars(parameters::Table& p, const std::string& mod="")
    : BASE{p,mod},
      seed(p.addTitle("StochasticTrajectory",mod).add("seed",mod,"Random number generator seed",1001ul)),
      prngStream(p.add("prngStream",mod,"Random number generator independent stream ordinal",1ul))
    {}
};


#ifndef NDEBUG
#pragma GCC warning "TODO: implement nextStream optionally from previousState"
#endif // NDEBUG
template<typename BASE>
pcg64 streamOfOrdo(const Pars<pcg64,BASE>& p/*, std::optional<EngineState> previousState = std::nullopt*/)
{
  return pcg64{p.seed,1001+p.prngStream};
}


template<typename BASE>
void incrementForNextStream(const Pars<pcg64,BASE>& p) {++p.prngStream;}


} // randomutils


namespace boost { namespace serialization {


template<typename Ar>
void load(Ar& ar, XoshiroCpp::Xoshiro256PlusPlus& xpp, unsigned) {
  XoshiroCpp::Xoshiro256PlusPlus::state_type state;
  ar & state;
  xpp.deserialize(state);
}

template<typename Ar>
void save(Ar& ar, XoshiroCpp::Xoshiro256PlusPlus const& xpp, unsigned) {
  auto state{xpp.serialize()};
  ar & state;
}

template<typename Ar>
void serialize(Ar& ar, XoshiroCpp::Xoshiro256PlusPlus& xpp, unsigned version) {
  if (typename Ar::is_saving())
    save(ar, xpp, version);
  else
    load(ar, xpp, version);
}


#define CPPQEDCORE_UTILS_RANDOM_H_REENTRANT
#define CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE std::mt19937_64
#define CPPQEDCORE_UTILS_RANDOM_H_STRING "mersenne_twister_engine state"

#include "Random.h"

#define CPPQEDCORE_UTILS_RANDOM_H_REENTRANT
#define CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE pcg64
#define CPPQEDCORE_UTILS_RANDOM_H_STRING "pcg64 state"

#include "Random.h"
  
} } // boost::serialization

#endif // CPPQEDCORE_UTILS_RANDOM_H_INCLUDED


#else // CPPQEDCORE_UTILS_RANDOM_H_REENTRANT


template<typename Ar>
void load(Ar& ar, CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE& randomEngine, unsigned) {
  std::string text;
  ar & text;
  std::istringstream iss(text);

  if (!(iss >> randomEngine))
    throw std::invalid_argument(CPPQEDCORE_UTILS_RANDOM_H_STRING);
}

template<typename Ar>
void save(Ar& ar, CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE const& randomEngine, unsigned) {
  std::ostringstream oss;
  if (!(oss << randomEngine))
    throw std::invalid_argument(CPPQEDCORE_UTILS_RANDOM_H_STRING);
  std::string text = oss.str();
  ar & text;
}

template<typename Ar>
void serialize(Ar& ar, CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE& randomEngine, unsigned version) {
  if (typename Ar::is_saving())
    save(ar, randomEngine, version);
  else
    load(ar, randomEngine, version);
}

#undef CPPQEDCORE_UTILS_RANDOM_H_STRING
#undef CPPQEDCORE_UTILS_RANDOM_H_RANDOMENGINE
#undef CPPQEDCORE_UTILS_RANDOM_H_REENTRANT

#endif // CPPQEDCORE_UTILS_RANDOM_H_REENTRANT
