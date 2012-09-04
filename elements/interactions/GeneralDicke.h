// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED

#include "GeneralDickeFwd.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode_.h"
#include "Spin.h"


class GeneralDicke : public structure::Interaction<2>, public structure::TridiagonalHamiltonian<2,true>
{
public:
  GeneralDicke(mode::SmartPtr, spin::SmartPtr, dcomp u, dcomp y);

};



#endif // ELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED
