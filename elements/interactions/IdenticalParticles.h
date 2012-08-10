// -*- C++ -*-
#ifndef   IDENTICAL_PARTICLES_INCLUDED
#define   IDENTICAL_PARTICLES_INCLUDED

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

  IdenticalParticles(const ParticleBase& part) : Base(ctorHelper(part)) {}

private:
  static const Particles ctorHelper(const ParticleBase& part) {Particles res; res=&part; return res;}

};


#endif // IDENTICAL_PARTICLES_INCLUDED
