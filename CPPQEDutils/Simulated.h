// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

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
template<typename ST, ode::system<ST> D, template <typename > typename OE> requires ode::engine<OE<std::decay_t<ST>>,std::decay_t<ST>>
struct Simulated
{
  Simulated(auto&& stateInit, D d, std::initializer_list<std::string> kl, OE<std::decay_t<ST>> oe)
    : state{std::forward<decltype(stateInit)>(stateInit)},
      derivs{d},
      keyLabels{kl},
      ode{oe} {
        if (keyLabels.size() < size(state)) keyLabels.insert(keyLabels.cend(),size(state)-keyLabels.size(),"N/A");
        else if (keyLabels.size() > size(state)) throw std::runtime_error("More keys than values in Simulated");
      }
  
  double time=0.;
  ST state;
  D derivs;
  std::list<std::string> keyLabels;
  OE<std::decay_t<ST>> ode;

  friend double getDtDid(const Simulated& s) {return getDtDid(s.ode);}

  friend double getTime(const Simulated& s) {return s.time;}

  friend LogTree logIntro(const Simulated& s) {return {{"Simulated",{"odeEngine",logIntro(s.ode)}}};}

  friend LogTree logOutro(const Simulated& s) {return logOutro(s.ode);}

  friend LogTree step(Simulated& s, double deltaT) {return step(s.ode,deltaT,s.derivs,s.time,s.state);}

  friend LogTree dataStreamKey(const Simulated& s) {return {{"Simulated",s.keyLabels}};}

  friend auto temporalDataPoint(const Simulated& s)
  {
    if constexpr (temporal_data_point<ST>) return s.state;
    else {
      return std::apply(hana::make_tuple,s.state);
    }
  }

  friend iarchive& readFromArrayOnlyArchive(Simulated& s, iarchive& iar) {return iar & s.state;}

  /** structure of Simulated archives:
  * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry )
  * state should precede time in order to be compatible with array-only archives
  */
  template <typename Archive>
  friend Archive& stateIO(Simulated& s, Archive& ar) {return stateIO(s.ode, ar & s.state & s.time);}

};


template <typename ST, typename D, template <typename > typename OE >
struct trajectory::MakeSerializationMetadata<Simulated<ST,D,OE>>
{
  /// This is probably not quite correct because elsewhere the typeID refers to the whole ST
  static auto _() {return SerializationMetadata{ boost::core::demangle( typeid(ST).name() ),
                                                 "Simulated",1};}
};


namespace simulated {

template<template <typename > typename OE, typename ST, typename D, typename ... ODE_EngineCtorParams>
requires ::cppqedutils::ode::engine<OE<std::decay_t<ST>>,std::decay_t<ST>>
auto make(ST&& stateInit, D derivs, std::initializer_list<std::string> keyLabels, ODE_EngineCtorParams&&... odePack)
{
  return Simulated<ST,D,OE>{std::forward<ST>(stateInit),derivs,keyLabels,OE<std::decay_t<ST>>{std::forward<ODE_EngineCtorParams>(odePack)...}};
}

template <typename BASE = Empty>
using Pars=::trajectory::Pars<::ode::Pars<BASE>>;


auto makeBoost(auto && stateInit, auto derivs, std::initializer_list<std::string> keyLabels, double dtInit, const auto& p)
{
  return make<ODE_EngineBoost>(std::forward<decltype(stateInit)>(stateInit),derivs,keyLabels,dtInit,p.epsRel,p.epsAbs);
}

} // simulated

