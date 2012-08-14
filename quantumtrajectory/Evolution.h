// -*- C++ -*-
#ifndef _EVOLUTION_H
#define _EVOLUTION_H

#include "impl/Composite.tcc"
#include "impl/EnsembleMCWF.tcc"
#include "impl/Evolution.tcc"
#include "ParsEvolution.h"
#include "impl/Master.tcc"
#include "impl/MCWF_Trajectory.tcc"
#include "QM_Picture.h"
#include "impl/StateVector.tcc"
#include "impl/Tridiagonal.tcc"

#include "BlitzArrayTraits.h"
#include "impl/EvolvedGSL.tcc"
#include "impl/Pars.tcc"

using parameters::ParameterTable    ;
using parameters::update            ;
using parameters::ParsNamedException;

#endif // _EVOLUTION_H
