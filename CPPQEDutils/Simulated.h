// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Convenience header for straightforward usage of the trajectory::Simulated class for classical simulations}
#ifndef CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED
#define CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED

#include "BlitzArray.h"
#include "Simulated_.h"

#include "Pars.h"

using ParameterTable=parameters::Table;
using parameters::update;

using namespace cppqedutils;

using trajectory::Pars;


template <int RANK, typename Derivs, typename ODE_Engine>
struct trajectory::MakeSerializationMetadata<Simulated<DArray<RANK>,Derivs,ODE_Engine>>
{
  static auto _() {return SerializationMetadata{"DArray","Simulated",RANK};}
};


template <int RANK, typename Derivs, typename ODE_Engine>
struct trajectory::MakeSerializationMetadata<Simulated<CArray<RANK>,Derivs,ODE_Engine>>
{
  static auto _() {return SerializationMetadata{"CArray","Simulated",RANK};}
};


#endif // CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED
