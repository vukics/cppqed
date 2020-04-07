// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_IDENTICALPARTICLES_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_IDENTICALPARTICLES_H_INCLUDED

#include "IdenticalParticlesFwd.h"

#include "Particle_.h"

#include "Interaction.h"


template<int NP>
class IdenticalParticles : public structure::Interaction<NP>
// A dummy implementation
{
public:
  typedef structure::Interaction<NP> Base;
  typedef typename Base::Frees Particles;

  IdenticalParticles(particle::Ptr       part) : Base(part,part) {}

};


#endif // CPPQEDELEMENTS_INTERACTIONS_IDENTICALPARTICLES_H_INCLUDED
