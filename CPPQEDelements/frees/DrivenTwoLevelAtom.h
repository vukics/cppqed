// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

// This is a very primitive implementation used as a test of the
// framework: Structure, MCWF_Trajectory, Master, EnsembleMCWF,
// etc. If there is an error, it is probably in there and not
// here. The corresponding driver is DTLA_C++QED. Better
// implementation based on class composition is found in Qbit_.h

#pragma once

#include "Qbit_.h"

#include "ElementLiouvillian.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "Hamiltonian.h"

#include "CMatrix.h"




class DrivenTwoLevelAtom : public structure::Free, public structure::ElementLiouvillian<1,1>, public qbit::Averaged
{
public:
  const dcomp getZa() const {return za_;}

protected:
  DrivenTwoLevelAtom(const qbit::ParsDrivenDissipative&);

private:
  double rate(structure::NoTime, const qbit::LazyDensityOperator&) const;

  void doActWithJ(structure::NoTime, qbit::StateVectorLow&) const;

  const dcomp za_, eta_;

};


// In Schroedinger picture the Hamiltonian is implemented as a CMatrix
class DrivenTwoLevelAtomSch : public DrivenTwoLevelAtom, public structure::HamiltonianTimeDependenceDispatched<1,structure::TimeDependence::NO>
{
public:
  DrivenTwoLevelAtomSch(const qbit::ParsDrivenDissipative&);
  
private:
  void addContribution_v(structure::NoTime, const qbit::StateVectorLow& psi, qbit::StateVectorLow& dpsidt) const {linalg::apply(psi,dpsidt,hamiltonianOverI_);}

  static const linalg::CMatrix hamiltonianOverI(dcomp za, dcomp etat);

  const linalg::CMatrix hamiltonianOverI_;

};

