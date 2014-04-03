// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParticleCavity.h"
#include "ParticleCavity_InterferenceWorkaround.h"

//#include "ParsParticleCavity.h"
#include "MathExtensions.h"
#include "Mode.h"

#include<boost/assign/list_of.hpp>


using namespace std;
using namespace boost::assign;
using namespace mathutils;

using particle:: mfNKX;
using particle::cosNKX;

using namespace mode;

namespace particlecavity_interferenceworkaround {

InterferenceBase::InterferenceBase(mode::Ptr mode, const particle::Ptr particle, double u, size_t kCav, ModeFunctionType modeCav)
  : MF_Base(modeCav,kCav),
    structure::Interaction<2>(Frees(mode,particle),
			      CF{"u",u,mode->getDimension()}),
    TridiagonalHamiltonian(particlecavity::interferic(mode,particle,sqr(u),u,MF_Base::member))
{
  getParsStream()<<"# Interference term with "<<getMF()<<endl;
}

} // particlecavity_interferenceworkaround
