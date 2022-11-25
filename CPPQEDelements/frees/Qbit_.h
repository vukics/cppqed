// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the qbit-bundle (tackling the dynamics of a single qbit)}
#pragma once

#include "ParsQbit.h"

#include "QuantumSystemDynamics.h"


namespace qbit {

using namespace ::structure; using namespace ::quantumdata;

TimeIndependentTerm<1> diagonalH(dcomp z);
TimeIndependentTerm<1> offDiagonalH(dcomp eta);

UnaryDiagonalPropagator<> propagator(dcomp z);

TimeIndependentJump<1> sigmaJump(double gamma_m);
TimeIndependentJump<1> sigmaPlusJump(double gamma_p);
TimeIndependentJump<1> sigma_zJump(double gamma_phi);

TimeIndependentSuperoperator<1> sigmaSuperoperator(double gamma_m);
TimeIndependentSuperoperator<1> sigmaPlusSuperoperator(double gamma_p);
TimeIndependentSuperoperator<1> sigma_zSuperoperator(double gamma_p);

TimeIndependentRate<1> sigmaRate(double gamma_m);
TimeIndependentRate<1> sigmaPlusRate(double gamma_p);
TimeIndependentRate<1> sigma_zRate(double gamma_p);

Lindblad<1> loss(double gamma_m) {return {"loss", sigmaJump(gamma_m), sigmaRate(gamma_m), sigmaSuperoperator(gamma_m)};}
Lindblad<1> gain(double gamma_p) {return {"gain", sigmaPlusJump(gamma_p), sigmaPlusRate(gamma_p), sigmaPlusSuperoperator(gamma_p)};}
Lindblad<1> dephasing(double gamma_phi) {return {"dephasing", sigma_zJump(gamma_phi), sigma_zRate(gamma_phi), sigma_zSuperoperator(gamma_phi)};}

static const ExpectationValue<1>
  population{.labels={"population"},
             .eva=[] (LazyDensityOperator<1> rho) {EV_Array res(1); res[0]=real(rho(1,1)); return res;}},
  polarization{.labels={"polarization"},
               .eva=[] (LazyDensityOperator<1> rho) {EV_Array res(2); dcomp pol{rho(1,0)}; res[0]=real(pol); res[1]=imag(pol); return res;}};


StateVector<1> state0();// {return mode::fock(0,2);}
StateVector<1> state1();// {return mode::fock(1,2);}
StateVector<1> init(dcomp psi1);


QuantumSystemDynamics<1> make(dcomp zSch, dcomp zI, dcomp eta, double gamma_m, double gamma_p, double gamma_phi)
{
  
}


} // qbit

