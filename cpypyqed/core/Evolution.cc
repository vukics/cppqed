// -*- C++ -*-

#include "PythonExtension.h"
#include "ParsPropertyMacros.h"

#include "BlitzArrayTraits.h"

#include "EvolvedGSL.tcc"
#include "Pars.h"

#include "QuantumSystem.h"

#include "StateVector.tcc"

#include "EnsembleMCWF.tcc"
#include "Evolution.tcc"
#include "Master.tcc"
#include "ParsEvolution.h"
#include "ParsMCWF_Trajectory.h"

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

using namespace boost::python;

using quantumtrajectory::ParsMCWF;
using trajectory::ParsRun;
using quantumdata::StateVector;
using structure::QuantumSystem;

namespace pythonext {

PARS_GETTER_SETTER(EvolutionMode, ParsEvolution, evol)
PARS_GETTER_SETTER(bool, ParsEvolution, negativity)



void export_25_Evolution()
{
  enum_<EvolutionMode>("EvolutionMode")
    .value("EM_SINGLE",      EM_SINGLE)
    .value("EM_ENSEMBLE",    EM_ENSEMBLE)
    .value("EM_MASTER",      EM_MASTER)
    .value("EM_MASTER_FAST", EM_MASTER_FAST)
  ;
  class_<ParsEvolution, bases<ParsRun,ParsMCWF> >
    (
      "ParsEvolution",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(evol)
    .PARS_PROPERTY(negativity)
  ;


#define EVOLVE_INSTANTIATIONS(z,r,data) \
  def("evolve", (void (*)(StateVector<r>&, const QuantumSystem<r>&, const ParsEvolution&)) &evolve);
BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_ADD(PYTHON_HALF_RANK,1), EVOLVE_INSTANTIATIONS, data)

}


} // pythonext