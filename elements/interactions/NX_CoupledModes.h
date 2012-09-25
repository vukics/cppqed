// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED

#include "NX_CoupledModesFwd.h"

#include "Mode_.h"

#include "Interaction.h"

namespace nxcoupledmodes {

class Base: public structure::Interaction<2>, public structure::TridiagonalHamiltonian<2,true>
{
public:
  Base(mode::Ptr, mode::Ptr, double u);

};

} // nxcoupledmodes



#define BIG_NAMESPACE_NAME             nxcoupledmodes
#define BIG_CLASS_NAME                 NX_CoupledModes
#define BIG_ADDITIONAL_PARAMETERS      , double u
#define BIG_ADDITIONAL_PARAMETERS_PASS ,u

#include "details/BinaryInteractionGenerator.h"

#endif // ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED
