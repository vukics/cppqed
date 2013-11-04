// -*- C++ -*-
#ifndef   INTERACTIONS_PARSMLJC_H_INCLUDED
#define   INTERACTIONS_PARSMLJC_H_INCLUDED

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


#include<ParsMLJC.tcc>


#endif // INTERACTIONS_PARSMLJC_H_INCLUDED
