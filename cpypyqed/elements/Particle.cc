// Copyright Raimar Sandner 2012–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "Free.h"
#include "Particle_.h"
#include "ParsParticle.h"
#include "QM_PictureFwd.h"
#include "StateVector.tcc"

using namespace boost::python;

using particle::Pars; using particle::ParsPumped;
using structure::freesystem::StateVector;

namespace pythonext {

namespace {
// Boost spits some ugly warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
BOOST_PYTHON_FUNCTION_OVERLOADS(wp_overloads,particle::wavePacket,1,2)
BOOST_PYTHON_FUNCTION_OVERLOADS(ho_overloads,particle::hoState,1,2)
#pragma GCC diagnostic pop
const StateVector (*wp1)(const Pars&, bool) = &particle::wavePacket;
const StateVector (*wp2)(const ParsPumped&, bool) = &particle::wavePacket;
const StateVector (*ho1)(const Pars&, bool) = &particle::hoState;
const StateVector (*ho2)(const ParsPumped&, bool) = &particle::hoState;
}

void export_Particle()
{
  class_<ParticleBase, bases<structure::QuantumSystem<1>>, boost::noncopyable>("ParticleBase", no_init);
  class_<PumpedParticleBase, bases<ParticleBase>,boost::noncopyable>("PumpedParticleBase", no_init);

  {
    scope namespaceScope = particleNameSpace;
    def("make", &particle::make, with_custodian_and_ward_postcall<0,1>());
    particleNameSpace.staticmethod("make");
    def("makePumped", &particle::makePumped, with_custodian_and_ward_postcall<0,1>());
    particleNameSpace.staticmethod("makePumped");
    register_ptr_to_python< particle::Ptr >();
    register_ptr_to_python< particle::PtrPumped >();
    implicitly_convertible<boost::shared_ptr<ParticleBase>,particle::Ptr>();
    implicitly_convertible<boost::shared_ptr<PumpedParticleBase>,particle::PtrPumped>();

    def("wavePacket", wp1, wp_overloads());
    def("wavePacket", wp2, wp_overloads());
    particleNameSpace.staticmethod("wavePacket");
    def("hoState", ho1, ho_overloads());
    def("hoState", ho2, ho_overloads());
    particleNameSpace.staticmethod("hoState");
    def("init", particle::init);
    particleNameSpace.staticmethod("init");
  }
}

} // pythonext