// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "KeyPrinter.h"
#include "Trajectory.h"

#include <boost/core/demangle.hpp>


using namespace cppqedutils;

/// Class fully implementing the trajectory::adaptive and trajectory::archiving concepts by streaming (and serializing) the whole content of the evolved array
/**
 * Realizes a simple ODE evolution with intermittent streamings
 * 
 * <b>Example usage:</b> simulation of a complex driven damped harmonic oscillator mode described by the ODE \f[\ddot{y}+2\gamma\,\dot{y}+y=e^{i\,\omega t},\f]
 * where \f$\gamma\f$ is the damping rate and \f$\omega\f$ the driving frequency, and the timescale has been chosen such that the eigenfrequency is 1.
 * 
 * \include HarmonicOscillatorComplex.cc
 */
template<typename ST, ode::system<ST> D, ode::engine<ST,D> OE>
struct Simulated
{
  using StreamedArray=ST;
  
  Simulated(auto&& stateInit, D d, std::initializer_list<std::string> keyLabels, OE oe)
    : state{std::forward<decltype(stateInit)>(stateInit)},
      derivs{d},
      keyPrinter{"Simulated",keyLabels},
      ode{oe} {
        auto& labels{keyPrinter.getLabels()};
        if (labels.size() < size(state)) labels.insert(labels.end(),size(state)-labels.size(),"N/A");
        else if (labels.size() > size(state)) throw std::runtime_error("More keys than values in Simulated");
      }
  
  double time=0.;
  ST state;
  D derivs;
  KeyPrinter keyPrinter;
  OE ode;
  
};


template<typename ST, typename D, typename OE>
double getTime(const Simulated<ST,D,OE>& s) {return s.time;}

template<typename ST, typename D, typename OE>
double getDtDid(const Simulated<ST,D,OE>& s) {return getDtDid(s.ode);}

template<typename ST, typename D, typename OE>
std::ostream& streamIntro(const Simulated<ST,D,OE>& s, std::ostream& os) {return streamIntro(s.ode,os<<"\nSimulated.\n");}

template<typename ST, typename D, typename OE>
std::ostream& streamOutro(const Simulated<ST,D,OE>& s, std::ostream& os) {return streamOutro(s.ode,os);}

template<typename ST, typename D, typename OE>
void step(Simulated<ST,D,OE>& s, double deltaT, std::ostream& logStream) {step(s.ode,deltaT,logStream,s.derivs,s.time,s.state);}

template<typename ST, typename D, typename OE>
std::ostream& streamKey(const Simulated<ST,D,OE>& s, std::ostream& os) {size_t i=3; return s.keyPrinter.stream(os,i);}

template<typename ST, typename D, typename OE>
ST stream(const Simulated<ST,D,OE>& s, std::ostream& os, int precision)
{
  for (size_t i=0; i < size(s.state); i++) os<<FormDouble(precision)(s.state[i])<<' ';
  return s.state;
}

template<typename ST, typename D, typename OE>
iarchive& readFromArrayOnlyArchive(Simulated<ST,D,OE>& s, iarchive& iar) {return iar & s.state;}

/** structure of Simulated archives:
 * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry )
 * state should precede time in order to be compatible with array-only archives
 */
template <typename ST, typename D, typename OE, typename Archive>
Archive& stateIO(Simulated<ST,D,OE>& s, Archive& ar) {return stateIO(s.ode, ar & s.state & s.time);}



template <typename ST, typename D, typename OE>
struct trajectory::MakeSerializationMetadata<Simulated<ST,D,OE>>
{
  /// This is probably not quite correct because elsewhere the typeID refers to the whole ST
  static auto _() {return SerializationMetadata{ boost::core::demangle( typeid(ST).name() ),
                                                 "Simulated",1};}
};


namespace simulated {

template<typename OE, typename ST, typename D, typename ... ODE_EngineCtorParams>
auto make(ST&& stateInit, D derivs, std::initializer_list<std::string> keyLabels, ODE_EngineCtorParams&&... odePack)
{
  return Simulated<std::decay_t<ST>,D,OE>{std::forward<ST>(stateInit),derivs,keyLabels,OE{std::forward<ODE_EngineCtorParams>(odePack)...}};
}

template <typename BASE = Empty>
using Pars=::trajectory::Pars<::ode::Pars<BASE>>;


template<typename ST, typename D, typename P>
auto makeBoost(ST&& stateInit, D derivs, std::initializer_list<std::string> keyLabels, double dtInit, const P& p)
{
  return make<ODE_EngineBoost<std::decay_t<ST>>>(std::forward<ST>(stateInit),derivs,keyLabels,dtInit,p.logControl,p.epsRel,p.epsAbs);
}

} // simulated

