// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// Services for dealing with frequency-like parameters, real or complex, for all physical systems
/** 
 * Such parameters need a special treatment because every such parameter from all the subsystems (either frees or interactions) of the given physical system has to be considered
 * as a possible largest frequency of the whole system. This largest frequency is needed for determining the initial time-step of the ODE routine.
 * 
 * Moreover, these parameters (together with all other, non-frequency-like parameters) have to be communicated towards the user when the framework
 * summarizes the parameters of a given run. Therefore, for each such parameter, the class stores not only the value, but also the name of the parameter,
 * plus another real number which multiplies the value of the named frequency, to give the actual frequency as appearing in the ODE.
 * 
 * Example: Consider the Hamiltonian of a free mode: \f[\omega a^\dagger a.\f]
 * In this case, the parameter supplied by the user is \f$\omega\f$, but the largest frequency appearing in the ODE is actually this times the dimension of the system
 * (the cutoff of the Fock space). Hence, the tuple stored by the class for this particular frequency-like parameter will be:
 * ~~~
 * ("omega",omega,cutoff)
 * ~~~
 * 
 * In the case of a pumping term of the Hamiltonian: \f[\eta\lp a^\dagger+a\rp,\f] the multiplier will be different:
 * ~~~
 * "eta",eta,sqrt(cutoff)
 * ~~~
 */
#pragma once

#include "ExpectationValues.h"
#include "Liouvillian.h"

#include <list>


namespace structure {


using SystemFrequencyDescriptor = std::tuple<std::string,std::variant<dcomp,double>,double>; /// name-value-multiplier

using SystemFrequencyStore = std::list<SystemFrequencyDescriptor>;


/// TODO: what happens if sfs is empty?
inline double highestFrequency(const SystemFrequencyStore& sfs)
{
  return sfs.size()
  ? std::ranges::max(sfs | std::views::transform( [](const auto& p) {
    return std::visit([&] (auto v) {
      return std::abs(v*get<2>(p));},get<1>(p)); } ) )
  : 0.;
}


template <typename T, size_t RANK>
concept quantum_system_dynamics = requires (T&& qsd)
{
  { getFreqs(qsd) } -> std::convertible_to<SystemFrequencyStore>;
  { getHa(qsd) } -> hamiltonian<RANK>;
  { getEx(qsd) } -> exact_propagator<RANK>;
  { getEV(qsd) } -> expectation_values<RANK>;
  { getLi(qsd) } -> std::convertible_to<Liouvillian<RANK>>;
  // { streamParameters(qsd,os) } -> std::convertible_to<std::ostream&>;
};


/// a simple implementation of the concept
/**
 * It’s probably not necessary to have this class at all.
 * Systems like Qubit, Mode could simply define a class for themselves that fulfill the quantum_system_dynamics concept
 * Also, in the concept, the getHa, getLi, getEV functions should be optional
 */
template <size_t RANK, hamiltonian<RANK> HA, exact_propagator<RANK> EX, expectation_values<RANK> EV>
struct QuantumSystemDynamics
{
  SystemFrequencyStore freqs;

  Liouvillian<RANK> li;

  HA ha; EX ex; EV ev;

  friend const SystemFrequencyStore& getFreqs(const QuantumSystemDynamics& qsd) {return qsd.freqs;}

  friend const Liouvillian<RANK>& getLi(const QuantumSystemDynamics& qsd) {return qsd.li;}

  friend const HA & getHa(const QuantumSystemDynamics& qsd) {return qsd.ha;}
  friend const EX & getEx(const QuantumSystemDynamics& qsd) {return qsd.ex;}
  friend const EV & getEV(const QuantumSystemDynamics& qsd) {return qsd.ev;}

};


template <size_t RANK, typename HA, typename EX, typename EV>
QuantumSystemDynamics(const SystemFrequencyStore&, const Liouvillian<RANK>&, HA&&, EX&&, EV&&) -> QuantumSystemDynamics<RANK,HA,EX,EV>;



// this class could be derived from QuantumSystemDynamics, with the latter having some kind of base class?
template <
  quantum_system_dynamics<1> QSD0,
  quantum_system_dynamics<1> QSD1,
  hamiltonian<2> HA,
  exact_propagator<2> EX,
  expectation_values<2> EV >
struct BinarySystem
{
  QSD0 qsd0; QSD1 qsd1;

  HA ha; EX ex; EV ev; // these are the properties that the interaction element might have


  BinarySystem(auto&& qsd0, size_t dim0, auto&& qsd1, size_t dim1, const SystemFrequencyStore& freqs, const Liouvillian<2> li, auto&& ha, auto&& ex, auto&& ev)
    : qsd0{std::forward<decltype(qsd0)>(qsd0)},
      qsd1{std::forward<decltype(qsd1)>(qsd1)},
      ha{std::forward<decltype(ha)>(ha)},
      ex{std::forward<decltype(ex)>(ex)},
      ev{std::forward<decltype(ev)>(ev)},
      freqsFull_{std::ranges::views::join({getFreqs(qsd0),getFreqs(qsd1),freqs})},
      liFull_{ [&] {
        std::vector< Lindblad<2> > res;
      } () },
      offsets0_{calculateSlicesOffsets<rA0>({dim0,dim1})},
      offsets1_{calculateSlicesOffsets<rA1>({dim0,dim1})}
  {}

  friend const auto& getFreqs(const BinarySystem& bs) {return
    std::vector<SystemFrequencyDescriptor> res;
    for (const auto& i : getFreqs(bs.qsd0)) res.push_back(i);
    for (const auto& i : getFreqs(bs.qsd1)) res.push_back(i);
    for (const auto& i : bs.freqs) res.push_back(i);
    return res;
  }


  friend auto getHa(const BinarySystem& bs) {
    return [&] (double t, StateVectorConstView<2> psi, StateVectorView<2> dpsidt, double t0) {
      hamiltonian_ns::broadcast<BinarySystem::rA0>(getHa(bs.qsd0),t,psi,dpsidt,t0,bs.offsets0_);
      hamiltonian_ns::broadcast<BinarySystem::rA1>(getHa(bs.qsd1),t,psi,dpsidt,t0,bs.offsets1_);
      bs.ha(t,psi,dpsidt,t0);
    };
  }

  friend auto getEx(const BinarySystem& bs) {
    return [&] (double t, lazy_density_operator<2> auto psi, double t0) {
      return hana::tuple(
        exact_propagator_ns::broadcast<BinarySystem::rA0>(getEx(bs.qsd0),t,psi,bs.offsets0_),
        exact_propagator_ns::broadcast<BinarySystem::rA1>(getEx(bs.qsd1),t,psi,bs.offsets1_),
        bs.ex(t,psi));
    };
  }

  /*
  friend void postProcessor(std::invoke_result_t<getEx,BinarySystem> )
  constexpr auto postProcessor(decltype(expectationValues)) {
  return [] (std::invoke_result_t<decltype(expectationValues),StateVectorConstView<1>> & t) {
    using namespace hana::literals; t[1_c] -= sqr(t[0_c]);
  };
  }*/

  friend const auto& getLi(const BinarySystem& bs) {return bs.li_;}

private:
  static constexpr auto
    rA0=retainedAxes<0>, rA1=retainedAxes<1>;

  const SystemFrequencyStore freqsFull_;

  const Liouvillian<2> liFull_;

  const std::vector<size_t> offsets0_, offsets1_;
};



} // structure

