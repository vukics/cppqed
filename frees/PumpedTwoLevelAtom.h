// -*- C++ -*-

// This is a very primitive implementation used as a test of the
// framework: Structure, MCWF_Trajectory, Master, EnsembleMCWF,
// etc. If there is an error, it is probably in there and not
// here. The corresponding driver is PTLA_C++QED. Better
// implementation based on class composition is found in Qbit.h

#ifndef _PUMPED_TWO_LEVEL_ATOM_H
#define _PUMPED_TWO_LEVEL_ATOM_H

#include "PumpedTwoLevelAtomFwd.h"

#include "QbitFwd.h"

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "TimeIndependentHamiltonian.h"

#include "CMatrix.h"

/*
#include "Algorithm.h"
#include "Pars.h"


*/



class PumpedTwoLevelAtom : public structure::Free, public structure::ElementLiouvillean<1>, public structure::ElementAveraged<1> 
{
public:
  typedef structure::ElementAveraged<1> Averaged;

  typedef structure::ElementLiouvillean<1> Base;
  typedef Base::StateVectorLow StateVectorLow;
  typedef Base::LazyDensityOperator LazyDensityOperator;

  double probability(const LazyDensityOperator&) const;

  void doActWithJ(StateVectorLow&) const;

  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&) const;

  const dcomp getZa() const {return za_;}

protected:
  PumpedTwoLevelAtom(const qbit::ParsPumpedLossy&);

private:
  const dcomp za_, eta_;

};


// In Schroedinger picture the Hamiltonian is implemented as a CMatrix
class PumpedTwoLevelAtomSch : public PumpedTwoLevelAtom, public structure::TimeIndependentHamiltonian<1>
{
public:
  typedef PumpedTwoLevelAtom Base;
  typedef Base::StateVectorLow StateVectorLow;

  PumpedTwoLevelAtomSch(const qbit::ParsPumpedLossy&);
  
  void addContribution(const StateVectorLow& psi, StateVectorLow& dpsidt) const 
  {linalg::apply(psi,dpsidt,hamiltonianOverI_);}

private:
  static const linalg::CMatrix hamiltonianOverI(const dcomp& za, const dcomp& etat);

  const linalg::CMatrix hamiltonianOverI_;

};

#endif // _PUMPED_TWO_LEVEL_ATOM_H
