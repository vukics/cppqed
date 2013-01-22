// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_IDENTICALPARTICLES_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_IDENTICALPARTICLES_H_INCLUDED

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

  IdenticalParticles(const ParticleBase& part) : Base(ctorHelper(cpputils::sharedPointerize(part))) {}
  IdenticalParticles(particle::Ptr       part) : Base(ctorHelper                           (part) ) {}

private:
  static const Particles ctorHelper(particle::Ptr part) {Particles res; res=part; return res;}

};


#endif // ELEMENTS_INTERACTIONS_IDENTICALPARTICLES_H_INCLUDED