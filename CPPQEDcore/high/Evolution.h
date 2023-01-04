// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Convenience header collecting to one place most of the components of C++QED core needed for writing scripts (basically, the core part of the Level-1 interface)}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED

#include "Evolution_.h"
#include "QM_Picture.h"
#include "Tridiagonal.tcc"

#include "Pars.h"


/// Introduces ParameterTable into the global namespace to break ambiguity between update and parameters::update
using ParameterTable = parameters::Table;

using parameters::update, picture::updateWithPicture;


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED

