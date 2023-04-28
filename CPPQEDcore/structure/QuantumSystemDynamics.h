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


namespace structure {


using SystemFrequencyDescriptor = std::tuple<std::string,std::variant<dcomp,double>,double>; /// name-value-multiplier


template <typename T>
concept system_frequency_store = 
  std::ranges::forward_range<T> && 
  std::is_same_v<std::ranges::range_value_t<T>,SystemFrequencyDescriptor>;


/// TODO: what happens if sfs is empty?
template <system_frequency_store SFS>
double highestFrequency(const SFS& sfs)
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
  { getFreqs(qsd) } -> system_frequency_store;
  // or even
  // { GetFrequencies<T>{}(qsd) } -> system_frequency_store;
  // which gives even more flexibility
  { getHa(qsd) } -> hamiltonian<RANK>;
  { getLi(qsd) } -> liouvillian<RANK>;
  { getEV(qsd) } -> expectation_values<RANK>;
  // { streamParameters(qsd,os) } -> std::convertible_to<std::ostream&>;
};


/// a simple implementation of the concept:
template <size_t RANK, system_frequency_store SFS, liouvillian<RANK> LI, hamiltonian<RANK> HA, expectation_values<RANK> EV>
struct QuantumSystemDynamics
{
  SFS freqs; LI li; HA ha; EV ev;

  friend auto getFreqs(const QuantumSystemDynamics& qsd) {return qsd.freqs;}
  friend auto getLi   (const QuantumSystemDynamics& qsd) {return qsd.li;}
  friend auto getHa   (const QuantumSystemDynamics& qsd) {return qsd.ha;}
  friend auto getEV   (const QuantumSystemDynamics& qsd) {return qsd.ev;}
};


template <system_frequency_store SFS, typename LI, typename HA, typename EV>
QuantumSystemDynamics(SFS&& freqs, LI&& li, HA&& ha, EV&& ev) -> QuantumSystemDynamics<std::ranges::range_value_t<LI>::N_RANK,SFS,LI,HA,EV>;


} // structure

