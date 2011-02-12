// -*- C++ -*-
#ifndef   PARS_MULTI_LEVEL_JAYNES_CUMMINGS_INCLUDED
#define   PARS_MULTI_LEVEL_JAYNES_CUMMINGS_INCLUDED

#include "MLJCFwd.h"

#include "ParsFwd.h"


namespace mljc {

template<typename VC>
struct Pars
{
  VC& gs;

  Pars(parameters::ParameterTable&, const std::string& ="");

};

} // mljc


#include<impl/ParsMLJC.tcc>


#endif // PARS_MULTI_LEVEL_JAYNES_CUMMINGS_INCLUDED
