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

/// name-value-multiplier
using SystemFrequencyDescriptor = std::tuple<std::string,std::variant<dcomp,double>,double>;

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
  { getDimensions(qsd) } -> std::convertible_to<Dimensions<RANK>>;
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
  Dimensions<RANK> dim;

  SystemFrequencyStore freqs;

  Liouvillian<RANK> li;

  HA ha;
  EX ex = exact_propagator_ns::noOp;
  EV ev = expectation_values_ns::noOp;

  friend const SystemFrequencyStore& getFreqs(const QuantumSystemDynamics& qsd) {return qsd.freqs;}

  friend const Liouvillian<RANK>& getLi(const QuantumSystemDynamics& qsd) {return qsd.li;}

  friend const HA & getHa(const QuantumSystemDynamics& qsd) {return qsd.ha;}
  friend const EX & getEx(const QuantumSystemDynamics& qsd) {return qsd.ex;}
  friend const EV & getEV(const QuantumSystemDynamics& qsd) {return qsd.ev;}

  friend auto getDimensions(const QuantumSystemDynamics& qsd) {return qsd.dim;}

};


template <size_t RANK, typename HA, typename EX, typename EV>
QuantumSystemDynamics(Dimensions<RANK>, const SystemFrequencyStore&, const Liouvillian<RANK>&, HA&&, EX&&, EV&&) -> QuantumSystemDynamics<RANK,HA,EX,EV>;

template <typename HA, typename EX, typename EV>
QuantumSystemDynamics(size_t, const SystemFrequencyStore&, const Liouvillian<1>&, HA&&, EX&&, EV&&) -> QuantumSystemDynamics<1,HA,EX,EV>;


namespace binary {

auto expectation_values(const quantum_system_dynamics<1> auto& qsd0, const quantum_system_dynamics<1> auto& qsd1,
                        const std::vector<size_t>& offsets0, const std::vector<size_t>& offsets1, 
                        const expectation_values<2> auto& ev) {
  return [&] (double t, lazy_density_operator<2> auto psi) {
    return hana::make_tuple(
      partialTrace<retainedAxes<0>,2>(psi,offsets0,[&] (auto m) {return expectation_values_ns::calculate<1>(getEV(qsd0),t,m);},plusTDP{}),
      partialTrace<retainedAxes<1>,2>(psi,offsets1,[&] (auto m) {return expectation_values_ns::calculate<1>(getEV(qsd1),t,m);},plusTDP{}),
      expectation_values_ns::calculate<2>(ev,t,psi));
  };
}
  
} // binary


// this class could be derived from QuantumSystemDynamics, with the latter having some kind of base class?
template <
  quantum_system_dynamics<1> QSD0,
  quantum_system_dynamics<1> QSD1,
  hamiltonian<2> HA,
  exact_propagator<2> EX,
  expectation_values<2> EV >
class BinarySystem
{
private:
  const SystemFrequencyStore freqsFull_;

  const std::vector<size_t> offsets0_, offsets1_;

  const Liouvillian<2> liFull_;

public:
  QSD0 qsd0; QSD1 qsd1;

  HA ha; EX ex; EV ev; // these are the properties that the interaction element might have

  BinarySystem(auto&& qsd0, auto&& qsd1, const SystemFrequencyStore& freqs, const Liouvillian<2>& li, auto&& ha, auto&& ex, auto&& ev)
    : freqsFull_{ [&] {
        SystemFrequencyStore res(getFreqs(qsd0));/* res.append_range(getFreqs(qsd1)); res.append_range(freqs);*/ return res;
      } () },
      offsets0_{calculateSlicesOffsets<retainedAxes<0>>(concatenate(getDimensions(qsd0),getDimensions(qsd1)))},
      offsets1_{calculateSlicesOffsets<retainedAxes<1>>(concatenate(getDimensions(qsd0),getDimensions(qsd1)))},
      liFull_{ [&] {
        Liouvillian<2> res(size(getLi(qsd0))+size(getLi(qsd1))+size(li));
        auto resIter=res.begin();
        for (const Lindblad<1> & l : getLi(qsd0) ) *resIter++ = liouvillian_ns::broadcast<retainedAxes<0>,2>(l,offsets0_);
        for (const Lindblad<1> & l : getLi(qsd1) ) *resIter++ = liouvillian_ns::broadcast<retainedAxes<1>,2>(l,offsets1_);
        for (const Lindblad<2> & l : li ) *resIter++ = l;
        return res;
      } () },
      qsd0{std::forward<decltype(qsd0)>(qsd0)},
      qsd1{std::forward<decltype(qsd1)>(qsd1)},
      ha{std::forward<decltype(ha)>(ha)},
      ex{std::forward<decltype(ex)>(ex)},
      ev{std::forward<decltype(ev)>(ev)}
  {}


  friend const SystemFrequencyStore& getFreqs(const BinarySystem& bs) {return bs.freqsFull_;}

  friend auto getHa(const BinarySystem& bs) {
    return [&] (double t, StateVectorConstView<2> psi, StateVectorView<2> dpsidt, double t0) {
      hamiltonian_ns::broadcast<retainedAxes<0>>(getHa(bs.qsd0),t,psi,dpsidt,t0,bs.offsets0_);
      hamiltonian_ns::broadcast<retainedAxes<1>>(getHa(bs.qsd1),t,psi,dpsidt,t0,bs.offsets1_);
      applyHamiltonian(bs.ha,t,psi,dpsidt,t0);
    };
  }

  // friend LogTree label(decltype(
  //   binary::hamiltonian( std::declval<QSD0>(), std::declval<QSD1>(), std::vector<size_t>{}, std::vector<size_t>{}, std::declval<HA>() )
  // ) ) {return "BinarySystem";}

  friend auto getEx(const BinarySystem& bs) {return exact_propagator_ns::noOp;}

  friend auto getEV(const BinarySystem& bs) {return binary::expectation_values(bs.qsd0,bs.qsd1,bs.offsets0_,bs.offsets1_,bs.ev);}
  
  // friend void postProcessor( decltype( getEV( std::declval<BinarySystem>() ) ) ) ;
  
  
/*
  friend auto getEx(const BinarySystem& bs) {
    return [&] (double t, StateVectorView<2> psi, double t0) {
      exact_propagator_ns::broadcast<BinarySystem::retainedAxes<0>>(getEx(bs.qsd0),t,psi,bs.offsets0_);
      exact_propagator_ns::broadcast<BinarySystem::retainedAxes<1>>(getEx(bs.qsd1),t,psi,bs.offsets1_);
      applyPropagator(bs.ex,t,psi,t0);
    };
  }

  friend auto getEV(const BinarySystem& bs) {
    return [&] (double t, lazy_density_operator<2> auto psi) {
      return hana::make_tuple(
        partialTrace<BinarySystem::retainedAxes<0>,2>(psi,bs.offsets0_,[&] (lazy_density_operator<1> auto m) {return getEV(bs.qsd0)(t,m);}),
        partialTrace<BinarySystem::retainedAxes<1>,2>(psi,bs.offsets1_,[&] (lazy_density_operator<1> auto m) {return getEV(bs.qsd1)(t,m);}),
        bs.ev(t,psi));
    };
  }

  friend LogTree label(const decltype(getEV(std::declval<BinarySystem>())) & ) {return "BinarySystem";}

  friend void postProcessor(std::invoke_result_t<getEx,BinarySystem> )
  constexpr auto postProcessor(decltype(expectationValues)) {
  return [] (std::invoke_result_t<decltype(expectationValues),StateVectorConstView<1>> & t) {
    using namespace hana::literals; t[1_c] -= sqr(t[0_c]);
  };
  }*/

  friend const auto& getLi(const BinarySystem& bs) {return bs.liFull_;}

  friend auto getDimensions(const BinarySystem& bs) {return concatenate(getDimensions(bs.qsd0),getDimensions(bs.qsd1));}

};



template <typename BS >
LogTree label( std::invoke_result_t< decltype(BS::getHa), BS > ) {return "BinarySystem";}



template <typename QSD0, typename QSD1, typename HA, typename EX, typename EV>
BinarySystem(QSD0, QSD1, const SystemFrequencyStore&, const Liouvillian<2>&, HA&& ha, EX&& ex, EV&& ev) -> BinarySystem<QSD0,QSD1,HA,EX,EV>;


} // structure

