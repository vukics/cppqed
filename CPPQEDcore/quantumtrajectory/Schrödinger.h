// Copyright András Vukics 2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumSystemDynamics.h"
#include "QuantumTrajectory.h"


namespace quantumtrajectory {


namespace schrödinger { typedef ode::Pars<> Pars; } // schrödinger


template <size_t RANK,
          ::structure::quantum_system_dynamics<RANK> QSD,
          ::cppqedutils::ode::engine<StorageType> OE>
struct Schrödinger
{
  Schrödinger(auto&& qsd, auto&& psi, auto&& oe)
    : qsd{std::forward<decltype(qsd)>(qsd)}, psi{std::forward<decltype(psi)>(psi)}, oe{std::forward<decltype(oe)>(oe)} {}

  double time=0., time0=0.;
  QSD qsd;
  StateVector<RANK> psi;
  OE oe;

  friend double getDtDid(const Schrödinger& s) {return getDtDid(s.oe);}

  friend double getTime(const Schrödinger& s) {return s.time;}

  /// TODO: put here the system-specific things
  friend LogTree logIntro(const Schrödinger& s)
  {
    return {{"Schrödinger",{{"odeEngine",logIntro(s.oe)},{"System","TAKE FROM SYSTEM"}}}};
    // streamCharacteristics(s.qsd,streamParameters(s.qsd,streamIntro(s.oe,os)<<"Solving Schrödinger equation."<<endl<<endl ) )<<endl;
  }

  friend LogTree logOutro(const Schrödinger& s) {return logOutro(s.oe);}

  friend LogTree step(Schrödinger& s, double deltaT)
  {
    return step(s.oe, std::min(s.oe.dtTry,deltaT), [&] (const StorageType& psiRaw, StorageType& dpsidtRaw, double t)
    {
      applyHamiltonian(getHa(s.qsd),t,
                       StateVectorConstView<RANK>{s.psi.extents,s.psi.strides,0,psiRaw},
                       StateVectorView<RANK>{s.psi.extents,s.psi.strides,0,dpsidtRaw.begin(),std::ranges::fill(dpsidtRaw,0)},
                       s.time0);
    },s.time,s.psi.dataStorage());
  }

  friend LogTree dataStreamKey(const Schrödinger& s) {return {{"Schrödinger","TAKE FROM SYSTEM"}};}

  friend auto temporalDataPoint(const Schrödinger& s)
  {
    return calculateAndPostprocess<RANK>( getEV(s.qsd), s.time, LDO<StateVector,RANK>(s.psi) );
  }

  friend iarchive& readFromArrayOnlyArchive(Schrödinger& s, iarchive& iar) {return iar & s.psi;}

  /** structure of Schrödinger archives:
  * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry )
  * state should precede time in order to be compatible with array-only archives
  */
  template <typename Archive>
  friend auto& stateIO(Schrödinger& s, Archive& ar)
  {
    stateIO(s.oe, ar & s.psi & s.time);
    // The following is a no-op in the case of state output, since time0 is equated with time @ the end of each step, however, it is an essential step for state input!
    s.time0=s.time;
    return ar;
  }


};


template<typename QSD, typename SV, typename OE>
Schrödinger(QSD , SV , OE ) -> Schrödinger< multiArrayRank_v<SV>, QSD, OE >;


namespace schrödinger {


template <template<typename> class OE, typename QSD, typename SV>
auto make(QSD&& qsd, SV&& state, const Pars& p)
{
  constexpr size_t RANK=multiArrayRank_v<std::decay_t<SV>>;
  using ODE=OE<StorageType>;

  double iDt=initialTimeStep(getFreqs(qsd)); // precalculate, since qsd gets forwarded (i.e., potentially moved)

  return Schrödinger<RANK,QSD,ODE>{std::forward<QSD>(qsd),std::forward<SV>(state),ODE{iDt,p.epsRel,p.epsAbs}};
}


} // schrödinger


} // quantumtrajectory


template <size_t RANK, typename QSD, typename ODE_Engine>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::Schrödinger<RANK,QSD,ODE_Engine>>
{
  static auto _() {return SerializationMetadata{"CArray","Schrödinger",RANK};}
};

