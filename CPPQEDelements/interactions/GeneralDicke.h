// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode_.h"
#include "Spin.h"


namespace generaldicke {

class Base : public structure::Interaction<2>, public quantumoperator::TridiagonalHamiltonian<2,true>
{
public:
  Base(mode::Ptr, spin::Ptr, dcomp u, dcomp y);

};

} // generaldicke


#define BIG_NAMESPACE_NAME             generaldicke
#define BIG_CLASS_NAME                 GeneralDicke
#define BIG_ADDITIONAL_PARAMETERS      , dcomp u, dcomp y
#define BIG_ADDITIONAL_PARAMETERS_PASS , u, y

#include "details_BinaryInteractionGenerator.h"




#endif // CPPQEDELEMENTS_INTERACTIONS_GENERALDICKE_H_INCLUDED
