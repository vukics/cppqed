// Copyright András Vukics 2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumTrajectory.h"


namespace quantumtrajectory {


template <size_t RANK,
          ::structure::quantum_system_dynamics<RANK> QSD,
          ::cppqedutils::ode::engine<StorageType> OE>
struct Schrödinger
{
  Schrödinger(auto&& qsd, auto&& psi, auto&& oe)

  double time=0., time0=0.;
  QSD qsd;
  StateVector psi;

  friend double getDtDid(const Simulated& s) {return getDtDid(s.ode);}

  friend double getTime(const Simulated& s) {return s.time;}

  friend LogTree logIntro(const Simulated& s) {return LogTree{{"Simulated",{{"odeEngine",logIntro(s.ode)},{"system parameters",s.parameters}}}};}

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


} // quantumtrajectory


template <size_t RANK, typename QSD, typename ODE_Engine>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::Schrödinger<RANK,QSD,ODE_Engine>>
{
  static auto _() {return SerializationMetadata{"CArray","Schrödinger",RANK};}
};

