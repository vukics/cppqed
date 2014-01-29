/// \briefFile{Convenience header for straightforward usage of the trajectory::Simulated class for classical simulations}
// -*- C++ -*-
#ifndef UTILS_SIMULATED_H_INCLUDED
#define UTILS_SIMULATED_H_INCLUDED

#include "BlitzArrayTraits.h"
#include "EvolvedGSL.tcc"
#include "Simulated.tcc"

#include "Pars.tcc"

#include <boost/bind.hpp>

using parameters::ParameterTable;
using parameters::update        ;
using parameters::NamedException;

using trajectory::Simulated;
using trajectory::Pars     ;

#endif // UTILS_SIMULATED_H_INCLUDED
