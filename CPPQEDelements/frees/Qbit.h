// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the qbit-bundle (tackling the dynamics of a single qbit)}
#pragma once

// #include "ParsQbit.h"

#include "QuantumSystemDynamics.h"


namespace qbit {

using namespace ::structure; using namespace ::quantumdata;

auto diagonalH(dcomp z) { return [=] (StateVectorConstView<1> psi, StateVectorView<1> dpsidt) {dpsidt(1)-=z*psi(1);}; }

auto offDiagonalH(dcomp eta) { return [=] (StateVectorConstView<1> psi, StateVectorView<1> dpsidt) {dpsidt(0)+=-conj(eta)*psi(1); dpsidt(1)+=eta*psi(0);}; }

::cppqedutils::LogTree label(decltype(diagonalH(dcomp{}))) { return {"diagonalH"}; }
::cppqedutils::LogTree label(decltype(offDiagonalH(dcomp{}))) { return {"offDiagonalH"}; }


UnaryDiagonalPropagator<> propagator(dcomp z);

TimeIndependentJump<1> sigmaJump(double gamma_m);
TimeIndependentJump<1> sigmaPlusJump(double gamma_p);
TimeIndependentJump<1> sigma_zJump(double gamma_phi);

TimeIndependentSuperoperator<1> sigmaSuperoperator(double gamma_m);
TimeIndependentSuperoperator<1> sigmaPlusSuperoperator(double gamma_p);
TimeIndependentSuperoperator<1> sigma_zSuperoperator(double gamma_p);

TimeIndependentRate<1> sigmaRate(double gamma_m);
TimeIndependentRate<1> sigmaPlusRate(double gamma_p);
/*TimeIndependentRate<1> sigma_zRate(double gamma_p);*/

Lindblad<1> loss(double gamma_m) {return {"loss", sigmaJump(gamma_m), sigmaRate(gamma_m), sigmaSuperoperator(gamma_m)};}
Lindblad<1> gain(double gamma_p) {return {"gain", sigmaPlusJump(gamma_p), sigmaPlusRate(gamma_p), sigmaPlusSuperoperator(gamma_p)};}
// Lindblad<1> dephasing(double gamma_phi) {return {"dephasing", sigma_zJump(gamma_phi), sigma_zRate(gamma_phi), sigma_zSuperoperator(gamma_phi)};}


static constexpr auto expectationValues = [] (lazy_density_operator<1> auto rho) { return hana::make_tuple(_(rho,0),_(rho,0,1)); };

::cppqedutils::LogTree label(decltype(expectationValues)) { return {"population","polarization"}; }


StateVector<1> state0();// {return mode::fock(0,2);}
StateVector<1> state1();// {return mode::fock(1,2);}
StateVector<1> init(dcomp psi1);

/*
auto make(dcomp zSch, dcomp zI, dcomp eta, double gamma_m, double gamma_p, double gamma_phi)
{
  return QuantumSystemDynamics{
    {{"zSch",zSch,1},{"zI",zI,1},{"η",eta,1},{"γ_m",gamma_m,1},{"γ_p",gamma_p,1},{"γ_phi",gamma_phi,1}},
    {{"diagSch",diagonalH(zSch)},{"diagI",propagator(zI)},{"drive",offDiagonalH(eta)}},
    {loss(gamma_m),gain(gamma_p),dephasing(gamma_phi)},
    expectationValues
  };
}
*/

} // qbit

