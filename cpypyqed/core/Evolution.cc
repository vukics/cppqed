// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "ParsPropertyMacros.h"

#include "BlitzArrayTraits.h"

#include "EvolvedGSL.tcc"
#include "Pars.h"

#include "QuantumSystem.h"

#include "StateVector.tcc"

#include "EnsembleMCWF.tcc"
#include "Evolution.tcc"
#include "Master.tcc"
#include "ParsMCWF_Trajectory.h"

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

using namespace boost::python;
using namespace evolution;

using trajectory::ParsRun;
using quantumdata::StateVector;
using structure::QuantumSystem;

namespace pythonext {

PARS_GETTER_SETTER(Method, Pars, evol)
PARS_GETTER_SETTER(bool, Pars, negativity)



void export_25_Evolution()
{
  {
    scope namespaceScope = evolutionNameSpace;
    enum_<Method>("Method", "Wrapper of :core:`evolution::Method`")
      .value("SINGLE",      SINGLE)
      .value("ENSEMBLE",    ENSEMBLE)
      .value("MASTER",      MASTER)
      .value("MASTER_FAST", MASTER_FAST)
    ;
    class_<Pars, bases<ParsRun,quantumtrajectory::mcwf::Pars> >
      (
        "Pars", "Wrapper of :core:`evolution::Pars`",
      init<parameters::ParameterTable&, optional<const std::string&> >()
      [with_custodian_and_ward<1,2>()]
      )
      .PARS_PROPERTY(evol)
      .PARS_PROPERTY(negativity)
    ;
  }

#define EVOLVE_INSTANTIATIONS(z,r,data) \
  def("evolve", (void (*)(StateVector<r>&, const QuantumSystem<r>&, const Pars&)) &evolve, "Wrapper of :core:`evolve` with RANK="#r);
BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_ADD(PYTHON_HALF_RANK,1), EVOLVE_INSTANTIATIONS, data)

}


} // pythonext