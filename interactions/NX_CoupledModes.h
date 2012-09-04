// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED

#include "NX_CoupledModesFwd.h"

#include "Mode_.h"

#include "Interaction.h"


class NX_CoupledModes: public structure::Interaction<2>, public structure::TridiagonalHamiltonian<2,true>
{
public:
  NX_CoupledModes(mode::SmartPtr, mode::SmartPtr, double u);

};


#endif // ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED
