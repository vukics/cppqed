// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "StateVector.h"


namespace structure {

template <typename H, size_t RANK>
concept time_independent_term = requires(H&& h, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt) { h(psi,dpsidt); };

template <typename H, size_t RANK>
concept one_time_dependent_term = requires(H&& h, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt) { h(t,psi,dpsidt); };

template <typename H, size_t RANK>
concept two_time_dependent_term = requires(H&& h, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0) { h(t,psi,dpsidt,t0); };

template <typename H, size_t RANK>
concept one_time_dependent_propagator = requires(H&& h, double t, StateVectorView<RANK> psi) { h(t,psi); };

template <typename H, size_t RANK>
concept two_time_dependent_propagator = requires(H&& h, double t, StateVectorView<RANK> psi, double t0) { h(t,psi,t0); };

template <typename H, size_t RANK>
concept term_or_propagator = 
  time_independent_term<H,RANK> || one_time_dependent_term<H,RANK> || two_time_dependent_term<H,RANK> || 
  one_time_dependent_propagator<H,RANK> || two_time_dependent_propagator<H,RANK> ;



/// a label and a term
/** the `label` becomes a `prefix:label` when the system becomes part of a more complex system */
template <typename H, size_t RANK>
concept hamiltonian_term = requires (H&& h) {
  { label(h) } -> std::convertible_to<std::string>;
  { termOrPropagator(h) } -> term_or_propagator<RANK>;
};


template <typename H, size_t RANK>
concept hamiltonian = hana_sequence<H> && !!hana::all_of(
  decltype(hana::transform(std::declval<H>(), hana::typeid_)){},
  []<class T>(T) { return hamiltonian_term<typename T::type,RANK>; });


/// applying a Hamiltonian is by default interpreted as |dpsidt>+=H|psi>
/** However, if indexing is done carefully, `psi` and `dpsidt` can refer to the same underlying data */
template <size_t RANK>
void applyHamiltonian(const hamiltonian<RANK> auto& ha, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  hana::for_each(hana::transform(ha, [] (auto h) {return func(h);} ), [&]<typename H> (H h) {
    if constexpr (time_independent_term  <H,RANK>) h(psi,dpsidt);
    else if      (one_time_dependent_term<H,RANK>) h(t-t0,psi,dpsidt);
    else if      (two_time_dependent_term<H,RANK>) h(t,psi,dpsidt,t0);
  });
}


/// applying a Propagator is interpreted as replacing |psi> with U|psi>
template <size_t RANK>
void applyPropagator(const hamiltonian<RANK> auto& ha, double t, StateVectorView<RANK> psi, double t0)
{
  hana::for_each(hana::transform(ha, [] (auto h) {return func(h);} ), [&]<typename H> (H h) {
    if constexpr (one_time_dependent_propagator<H,RANK>) h(t-t0,psi);
    else if      (two_time_dependent_propagator<H,RANK>) h(t,psi,t0);
  });
}



/// A unary propagator that assumes that the operator that transforms between the pictures is diagonal
template<bool IS_TWO_TIME=false>
struct UnaryDiagonalPropagator
{
  using Diagonal = std::valarray<dcomp>;

  using UpdateFunctional = std::conditional_t<IS_TWO_TIME,std::function<void(double,Diagonal&,double)>,std::function<void(double,Diagonal&)>>;
  
  UnaryDiagonalPropagator(size_t dim, UpdateFunctional updateDiagonal) : diagonal{dim}, updateDiagonal_{updateDiagonal} {}
  
  void operator()(double t, StateVectorView<1> psi, double t0) const requires (IS_TWO_TIME);// {if (t!=t_ || t0!=t0_) {updateDiagonal_(t_=t,diagonal,t0_=t0);} psi*=diagonal;}

  void operator()(double t, StateVectorView<1> psi) const requires (!IS_TWO_TIME);// {if (t!=t_) {updateDiagonal_(t_=t,diagonal);} psi*=diagonal;}

  mutable Diagonal diagonal;

private:
  mutable double t_=0, t0_=0;
  
  const UpdateFunctional updateDiagonal_;

};


template <size_t RANK, term_or_propagator<RANK> TOP>
struct HamiltonianTerm
{
  std::string label;
  TOP termOrPropagator;
  
  friend std::string label(HamiltonianTerm ht) {return ht.label;}
  friend TOP termOrPropagator(HamiltonianTerm ht) {return ht.termOrPropagator;}
};


} // structure

