// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "StateVector.h"

#include <valarray>


namespace structure {


/// TODO: this could be fused with ode::system, but let’s not go into that for the moment :)
template <typename H, typename StateIn, typename StateOut>
concept ode_derivative_time_independent_contribution = requires(H&& h, StateIn psi, StateOut dpsidt) { h(psi,dpsidt); };

template <typename H, typename StateIn, typename StateOut>
concept ode_derivative_time_dependent_contribution = requires(H&& h, double t, StateIn psi, StateOut dpsidt) { h(t,psi,dpsidt); };


template <typename StateIn, typename StateOut>
using ODE_derivativeTimeIndependentFunctional = std::function<void(StateIn psi, StateOut dpsidt)> ;

template <typename StateIn, typename StateOut>
using ODE_derivativeTimeDependentFunctional = std::function<void(double t, StateIn psi, StateOut dpsidt)> ;


template <size_t RANK> using TimeIndependentTerm = ODE_derivativeTimeIndependentFunctional<StateVectorConstView<RANK>,StateVectorView<RANK>>;
template <size_t RANK> using TimeDependentTerm = ODE_derivativeTimeDependentFunctional<StateVectorConstView<RANK>,StateVectorView<RANK>>;


namespace hamiltonian_ns {

template <typename H, size_t RANK>
concept time_independent_term = ode_derivative_time_independent_contribution<H,StateVectorConstView<RANK>,StateVectorView<RANK>>;

template <typename H, size_t RANK>
concept one_time_dependent_term = ode_derivative_time_dependent_contribution<H,StateVectorConstView<RANK>,StateVectorView<RANK>>;

template <typename H, size_t RANK>
concept two_time_dependent_term = requires (H&& h, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0) { h(t,psi,dpsidt,t0); };

template <typename H, size_t RANK>
concept term = time_independent_term<H,RANK> || one_time_dependent_term<H,RANK> || two_time_dependent_term<H,RANK>;

template <typename H, size_t RANK>
concept one_time_dependent_propagator = requires (H&& h, double t, StateVectorView<RANK> psi) { h(t,psi); };

template <typename H, size_t RANK>
concept two_time_dependent_propagator = requires (H&& h, double t, StateVectorView<RANK> psi, double t0) { h(t,psi,t0); };

template <typename H, size_t RANK>
concept propagator = one_time_dependent_propagator<H,RANK> || two_time_dependent_propagator<H,RANK> ;

template <typename H, size_t RANK> concept functional = term<H,RANK> || propagator<H,RANK>;

} // hamiltonian_ns



/// applying a Hamiltonian term is by default interpreted as |dpsidt>+=H|psi>/(i*hbar)
/** However, if indexing is done carefully, `psi` and `dpsidt` can refer to the same underlying data */
template <size_t RANK, hamiltonian_ns::functional<RANK> T>
void applyHamiltonian(const T& h, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  if constexpr      (hamiltonian_ns::time_independent_term  <T,RANK>) h(psi,dpsidt);
  else if constexpr (hamiltonian_ns::one_time_dependent_term<T,RANK>) h(t-t0,psi,dpsidt);
  else if constexpr (hamiltonian_ns::two_time_dependent_term<T,RANK>) h(t,psi,dpsidt,t0);
}


/// applying the propagators (interpreted as replacing |psi> with U|psi>)
template <size_t RANK, hamiltonian_ns::functional<RANK> T>
void applyPropagator(const T& h, double t, StateVectorView<RANK> psi, double t0)
{
  if constexpr      (hamiltonian_ns::one_time_dependent_propagator<T,RANK>) h(t-t0,psi);
  else if constexpr (hamiltonian_ns::two_time_dependent_propagator<T,RANK>) h(t,psi,t0);
}


template <typename H, size_t RANK>
concept hamiltonian = ::cppqedutils::labelled<H> && hamiltonian_ns::functional<H,RANK>;



/// A unary propagator that assumes that the operator that transforms between the pictures is diagonal
template<bool IS_TWO_TIME=false>
struct UnaryDiagonalPropagator
{
  using Diagonal = std::valarray<dcomp>;

  using UpdateFunctional = std::conditional_t<IS_TWO_TIME,std::function<void(double,Diagonal&,double)>,std::function<void(double,Diagonal&)>>;
  
  UnaryDiagonalPropagator(std::string l, size_t dim, UpdateFunctional updateDiagonal) : label(l), diagonal(dim), updateDiagonal_{updateDiagonal} {}
  
  void operator()(double t, StateVectorView<1> psi, double t0) const requires (IS_TWO_TIME);// {if (t!=t_ || t0!=t0_) {updateDiagonal_(t_=t,diagonal,t0_=t0);} psi*=diagonal;}

  void operator()(double t, StateVectorView<1> psi) const requires (!IS_TWO_TIME);// {if (t!=t_) {updateDiagonal_(t_=t,diagonal);} psi*=diagonal;}

  const std::string label;

  mutable Diagonal diagonal;

private:
  mutable double t_=0, t0_=0;
  
  const UpdateFunctional updateDiagonal_;

};


template <size_t RANK, hamiltonian_ns::functional<RANK> F>
struct HamiltonianElement
{
  std::string label;

  F functional;
  
  void operator()(double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0) const requires hamiltonian_ns::term<F,RANK>
  {
    applyHamiltonian(functional,t,psi,dpsidt,t0);
  };

  void operator()(double t, StateVectorView<RANK> psi, double t0) const requires hamiltonian_ns::propagator<F,RANK>
  {
    applyPropagator(functional,t,psi,t0);
  };

};



template <size_t RANK, hamiltonian_ns::functional<RANK> F>
auto makeHamiltonianElement(std::string label, F&& functional)
{
  return HamiltonianElement<RANK,F>{.label{label},.functional{std::forward<F>(functional)}};
}


// itself a hamiltonian
template <size_t RANK, hamiltonian<RANK>... H>
struct HamiltonianCollection
{
  hana::tuple<H...> collection;

  void operator()(double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0) const
  {
    hana::for_each(collection, [&] (const auto& h) {applyHamiltonian(h,t,psi,dpsidt,t0);});
  }

  void operator()(double t, StateVectorView<RANK> psi, double t0) const
  {
    hana::for_each(collection, [&] (const auto& h) {applyPropagator(h,t,psi,t0);});
  }

  friend ::cppqedutils::LogTree label(HamiltonianCollection hc) {
    cppqedutils::json res;
    hana::for_each(hc.collection,[&] (const auto& h) {res.push_back(::cppqedutils::getLabel<::cppqedutils::LogTree>(h));});
    return res;
  }

};


template <size_t RANK, hamiltonian<RANK>... H>
auto makeHamiltonianCollection(H&&... h) { return HamiltonianCollection<RANK,H...>{.collection{std::forward<H>(h)...}}; }


static_assert(hamiltonian<UnaryDiagonalPropagator<false>,1>);
static_assert(hamiltonian<UnaryDiagonalPropagator<true>,1>);



/*
template <typename H, size_t RANK>
concept hamiltonian = hana_sequence<H> && !!hana::all_of(
  decltype(hana::transform(std::declval<H>(), hana::typeid_)){},
  []<class T>(T) { return hamiltonian_element<typename T::type,RANK>; });

*/


} // structure

