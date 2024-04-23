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

/// TODO: this and the previous could be put into a precompiled library
inline auto SFS_toJSON(const SystemFrequencyStore& sfs)
{
  LogTree res;
  for (const SystemFrequencyDescriptor& sfd : sfs)
    res.emplace(std::get<0>(sfd),std::visit(overload{
      [] (double v) {return json::value{v};},
      [] (dcomp v) {return json::value(json::array{v.real(),v.imag()});}
    },std::get<1>(sfd)));
  return res;
}


template <typename T, size_t RANK>
concept quantum_system_dynamics = labelled<T,std::string> && requires (T&& qsd)
{
  { getFreqs(qsd) } -> std::convertible_to<SystemFrequencyStore>;
  { getHa(qsd) } -> hamiltonian<RANK>;
  { getEx(qsd) } -> exact_propagator<RANK>;
  { getEV(qsd) } -> expectation_values<RANK>;
  { getLi(qsd) } -> std::convertible_to<Liouvillian<RANK>>;
  { getDimensions(qsd) } -> std::convertible_to<Dimensions<RANK>>;
  { getParameters(qsd) } -> std::convertible_to<LogTree>;
};


/// a simple implementation of the concept
/**
 * It’s probably not necessary to have this class at all.
 * Systems like Qubit, Mode could simply define a class for themselves that fulfill the quantum_system_dynamics concept
 * Also, in the concept, the getHa, getLi, getEV functions should be optional
 *
 * A JSON snippet could be supplied to express the other (non-frequency) parameters
 */
template <size_t RANK, hamiltonian<RANK> HA, exact_propagator<RANK> EX, expectation_values<RANK> EV>
struct QuantumSystemDynamics
{
  std::string label;

  Dimensions<RANK> dim;

  SystemFrequencyStore freqs;

  Liouvillian<RANK> li;

  HA ha;
  EX ex = exact_propagator_ns::noOp;
  EV ev = expectation_values_ns::noOp;

  LogTree nonFrequencyParams={};

  friend const SystemFrequencyStore& getFreqs(const QuantumSystemDynamics& qsd) {return qsd.freqs;}

  friend const Liouvillian<RANK>& getLi(const QuantumSystemDynamics& qsd) {return qsd.li;}

  friend const HA & getHa(const QuantumSystemDynamics& qsd) {return qsd.ha;}
  friend const EX & getEx(const QuantumSystemDynamics& qsd) {return qsd.ex;}
  friend const EV & getEV(const QuantumSystemDynamics& qsd) {return qsd.ev;}

  friend auto getDimensions(const QuantumSystemDynamics& qsd) {return qsd.dim;}

  friend auto getParameters(const QuantumSystemDynamics& qsd)
  {
    LogTree res{SFS_toJSON(getFreqs(qsd))};
    res.insert(qsd.nonFrequencyParams.cbegin(),qsd.nonFrequencyParams.cend());
    return res;
  }

};


template <size_t RANK, typename HA, typename EX, typename EV>
QuantumSystemDynamics(std::string, Dimensions<RANK>, const SystemFrequencyStore&, const Liouvillian<RANK>&, HA&&, EX&&, EV&&, const LogTree& = {})
-> QuantumSystemDynamics<RANK,HA,EX,EV>;

template <typename HA, typename EX, typename EV>
QuantumSystemDynamics(std::string, size_t, const SystemFrequencyStore&, const Liouvillian<1>&, HA&&, EX&&, EV&&, const LogTree& = {})
-> QuantumSystemDynamics<1,HA,EX,EV>;


namespace binary {

auto assembleEV(const quantum_system_dynamics<1> auto& qsd0, const quantum_system_dynamics<1> auto& qsd1,
                const std::vector<size_t>& offsets0, const std::vector<size_t>& offsets1,
                const expectation_values<2> auto& ev) {
  return [&] (double t, lazy_density_operator<2> auto psi) {
    return hana::make_tuple(
      partialTrace<retainedAxes<0>,2>(psi,offsets0,[&] (auto m) {return calculateExpectationValues<1>(getEV(qsd0),t,m);},plusTDP{}),
      partialTrace<retainedAxes<1>,2>(psi,offsets1,[&] (auto m) {return calculateExpectationValues<1>(getEV(qsd1),t,m);},plusTDP{}),
      calculateExpectationValues<2>(ev,t,psi));
  };
}


template <
  typename QSD0,
  typename QSD1,
  typename EV>
LogTree label(const decltype( assembleEV( std::declval<QSD0>(), std::declval<QSD1>(), std::vector<size_t>{}, std::vector<size_t>{}, std::declval<EV>() ) ) & ) {return {};}

} // binary



/**
 * TODO: this class should be templated for the full hamiltonian, propagator, expectation_values types
 * in that case, the ctor could calculate these functions, and BinarySystem expose them as members, instead of these getter-function craziness
 * (the quantum_system_dynamics has to be redefined for this)
 * the maker function / deduction guide can calculate these types rather easily
 */
template <
  quantum_system_dynamics<1> QSD0,
  quantum_system_dynamics<1> QSD1,
  hamiltonian<2> HA,
  exact_propagator<2> EX,
  expectation_values<2> EV,
  typename EVFULL>
class BinarySystem
{
private:
  const std::vector<size_t> offsets0_, offsets1_;

  const Liouvillian<2> liFull_;

  const EVFULL evFull_;

public:
  QSD0 qsd0; QSD1 qsd1;

  HA ha; EX ex; EV ev; // these are the properties that the interaction element might have

  const std::string label="BinarySystem";

  BinarySystem(auto&& qsd0, auto&& qsd1, const SystemFrequencyStore& freqs, const Liouvillian<2>& li, auto&& ha, auto&& ex, auto&& ev)
    : offsets0_{calculateSlicesOffsets<retainedAxes<0>>(concatenate(getDimensions(qsd0),getDimensions(qsd1)))},
      offsets1_{calculateSlicesOffsets<retainedAxes<1>>(concatenate(getDimensions(qsd0),getDimensions(qsd1)))},
      liFull_{ [&] {
        Liouvillian<2> res(size(getLi(qsd0))+size(getLi(qsd1))+size(li));
        auto resIter=res.begin();
        for (const Lindblad<1> & l : getLi(qsd0) ) *resIter++ = liouvillian_ns::broadcast<retainedAxes<0>,2>(l,offsets0_);
        for (const Lindblad<1> & l : getLi(qsd1) ) *resIter++ = liouvillian_ns::broadcast<retainedAxes<1>,2>(l,offsets1_);
        for (const Lindblad<2> & l : li ) *resIter++ = l;
        return res;
      } () },
      evFull_{binary::assembleEV(qsd0,qsd1,offsets0_,offsets1_,ev)},
      qsd0{std::forward<decltype(qsd0)>(qsd0)},
      qsd1{std::forward<decltype(qsd1)>(qsd1)},
      ha{std::forward<decltype(ha)>(ha)},
      ex{std::forward<decltype(ex)>(ex)},
      ev{std::forward<decltype(ev)>(ev)}
  {}


  friend auto getFreqs(const BinarySystem& bs) {
    SystemFrequencyStore res(getFreqs(bs.qsd0));
    { // TODO: append_range can be used in C++23
      const auto& f=getFreqs(bs.qsd1);
      res.insert(res.end(),f.begin(),f.end());
    }
//    res.append_range(freqs);
    return res;
  }

  friend auto getHa(const BinarySystem& bs) {
    return [&] (double t, StateVectorConstView<2> psi, StateVectorView<2> dpsidt, double t0) {
      hamiltonian_ns::broadcast<retainedAxes<0>>(getHa(bs.qsd0),t,psi,dpsidt,t0,bs.offsets0_);
      hamiltonian_ns::broadcast<retainedAxes<1>>(getHa(bs.qsd1),t,psi,dpsidt,t0,bs.offsets1_);
      applyHamiltonian(bs.ha,t,psi,dpsidt,t0);
    };
  }

  friend auto getEx(const BinarySystem& bs) {return exact_propagator_ns::noOp;}

  friend auto getEV(const BinarySystem& bs) {return bs.evFull_;}
  
  
/*
  friend auto getEx(const BinarySystem& bs) {
    return [&] (double t, StateVectorView<2> psi, double t0) {
      exact_propagator_ns::broadcast<BinarySystem::retainedAxes<0>>(getEx(bs.qsd0),t,psi,bs.offsets0_);
      exact_propagator_ns::broadcast<BinarySystem::retainedAxes<1>>(getEx(bs.qsd1),t,psi,bs.offsets1_);
      applyPropagator(bs.ex,t,psi,t0);
    };
  }


*/

//  friend LogTree label(const EVFULL&) {return {};}

  friend const auto& getLi(const BinarySystem& bs) {return bs.liFull_;}

  friend auto getDimensions(const BinarySystem& bs) {return concatenate(getDimensions(bs.qsd0),getDimensions(bs.qsd1));}

  friend LogTree getParameters(const BinarySystem& bs)
  {
    return {{getLabel(bs.qsd0),getParameters(bs.qsd0)},{getLabel(bs.qsd1),getParameters(bs.qsd1)}};
    // Interaction element could be just a tag class, but also have some functionalities, like checking Free subsystems for compatibility (e.g. Jaynes-Cummings expects a mode and a qbit)
  }
};



template <typename QSD0, typename QSD1, typename HA, typename EX, typename EV>
BinarySystem(QSD0&& qsd0, QSD1&& qsd1, const SystemFrequencyStore& sfs, const Liouvillian<2>& li, HA&& ha, EX&& ex, EV&& ev)
-> BinarySystem<QSD0,QSD1,HA,EX,EV,decltype(binary::assembleEV(qsd0,qsd1,std::vector<size_t>{},std::vector<size_t>{},ev))>;


} // structure

