// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_H_INCLUDED

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


#endif // CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_H_INCLUDED
