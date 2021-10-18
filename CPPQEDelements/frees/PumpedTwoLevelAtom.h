// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

// This is a very primitive implementation used as a test of the
// framework: Structure, MCWF_Trajectory, Master, EnsembleMCWF,
// etc. If there is an error, it is probably in there and not
// here. The corresponding driver is PTLA_C++QED. Better
// implementation based on class composition is found in Qbit_.h

#ifndef CPPQEDELEMENTS_FREES_PUMPEDTWOLEVELATOM_H_INCLUDED
#define CPPQEDELEMENTS_FREES_PUMPEDTWOLEVELATOM_H_INCLUDED

#include "Qbit_.h"

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "Hamiltonian.h"

#include "CMatrix.h"




class PumpedTwoLevelAtom : public structure::Free, public structure::ElementLiouvillean<1,1>, public qbit::Averaged
{
public:
  const dcomp getZa() const {return za_;}

protected:
  PumpedTwoLevelAtom(const qbit::ParsPumpedLossy&);

private:
  double rate(structure::NoTime, const LazyDensityOperator&) const;

  void doActWithJ(structure::NoTime, qbit::StateVectorLow&) const;

  const dcomp za_, eta_;

};


// In Schroedinger picture the Hamiltonian is implemented as a CMatrix
class PumpedTwoLevelAtomSch : public PumpedTwoLevelAtom, public structure::HamiltonianTimeDependenceDispatched<1,structure::TimeDependence::NO>
{
public:
  PumpedTwoLevelAtomSch(const qbit::ParsPumpedLossy&);
  
private:
  void addContribution_v(structure::NoTime, const qbit::StateVectorLow& psi, qbit::StateVectorLow& dpsidt) const {linalg::apply(psi,dpsidt,hamiltonianOverI_);}

  static const linalg::CMatrix hamiltonianOverI(dcomp za, dcomp etat);

  const linalg::CMatrix hamiltonianOverI_;

};

#endif // CPPQEDELEMENTS_FREES_PUMPEDTWOLEVELATOM_H_INCLUDED
