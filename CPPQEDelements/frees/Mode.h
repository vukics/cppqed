// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumSystemDynamics.h"

#include "MultiDiagonal.h"


namespace mode {

using namespace ::structure; using namespace ::quantumdata;
using MultiDiagonal = ::quantumoperator::MultiDiagonal<1> ;


MultiDiagonal aOp(size_t cutoff);

MultiDiagonal aDagOp(size_t cutoff);

MultiDiagonal nOp(size_t cutoff);


double photonnumber(lazy_density_operator<1> auto matrix)
{
  double res=0;
  for (size_t n=0; n<matrix.extents[0]; ++n) res+=n*_(matrix,n);
  return res;
}


auto hamiltonian(size_t cutoff, dcomp z, double omegaKerr, dcomp eta)
{
  std::string label{"Mode monolithic MultiDiagonal:"};
  MultiDiagonal ham, a{aOp(cutoff)}, aDag{aDagOp(cutoff)};

  if (abs(z)) {label+=" linear term;"; ham-=z*nOp(cutoff);}
  if (omegaKerr) {label+=" Kerr term;"; ham+=(omegaKerr/1.i) * (aDag|aDag|a|a);}
  if (abs(eta)) {label+=" drive term;"; ham+=twoTimesImagPartOf(eta*aDag);}

  return makeHamiltonianElement<1>(label,std::move(ham));
    //-z*nOp(cutoff) + (omegaKerr/1.i) * (aDag|aDag|a|a) + twoTimesImagPartOf(eta*aDag)
}



TimeIndependentJump<1> aJump(double fact);
TimeIndependentJump<1> aDagJump(double fact);

TimeIndependentSuperoperator<1> aSuperoperator(double fact);
TimeIndependentSuperoperator<1> aDagSuperoperator(double fact);


Lindblad<1> photonLoss(double kappa, double nTh);
Lindblad<1> photonGain(double kappa, double nTh);


static constexpr auto expectationValues = [] (lazy_density_operator<1> auto rho) { return photonnumber(rho); };

::cppqedutils::LogTree label(decltype(expectationValues)) { return {"photon number"}; }



/*


UnaryDiagonalPropagator<> propagator(dcomp z);


TimeIndependentSuperoperator<1> sigmaSuperoperator(double gamma_m);
TimeIndependentSuperoperator<1> sigmaPlusSuperoperator(double gamma_p);
TimeIndependentSuperoperator<1> sigma_zSuperoperator(double gamma_p);

TimeIndependentRate<1> sigmaRate(double gamma_m);
TimeIndependentRate<1> sigmaPlusRate(double gamma_p);

Lindblad<1> loss(double gamma_m) {return {"loss", sigmaJump(gamma_m), sigmaRate(gamma_m), sigmaSuperoperator(gamma_m)};}
Lindblad<1> gain(double gamma_p) {return {"gain", sigmaPlusJump(gamma_p), sigmaPlusRate(gamma_p), sigmaPlusSuperoperator(gamma_p)};}
// Lindblad<1> dephasing(double gamma_phi) {return {"dephasing", sigma_zJump(gamma_phi), sigma_zRate(gamma_phi), sigma_zSuperoperator(gamma_phi)};}


static constexpr auto expectationValues = [] (lazy_density_operator<1> auto rho) { return hana::make_tuple(real(_(rho,0,0)),_(rho,0,1)); };

::cppqedutils::LogTree label(decltype(expectationValues)) { return {"population","polarization"}; }


StateVector<1> state0();// {return mode::fock(0,2);}
StateVector<1> state1();// {return mode::fock(1,2);}
StateVector<1> init(dcomp psi1);
*/

} // mode
