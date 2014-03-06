// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"
#include "blitz2numpy.tcc"
#include "utils.h"

#include "ParsPropertyMacros.h"

#include "BlitzArrayTraits.h"

#include "EvolvedGSL.tcc"
#include "Pars.h"

#include "QuantumSystem.h"

#include "DensityOperator.h"
#include "StateVector.tcc"
#include "LazyDensityOperator.tcc"

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

template<int RANK>
object py_evolve(quantumdata::StateVector<RANK>& psi,
              typename structure::QuantumSystem<RANK>::Ptr sys,
              const evolution::Pars& p)
{
  using namespace quantumdata;

  typedef StateVector        <RANK> SV;
  typedef DensityOperator    <RANK> DO;
  typedef LazyDensityOperator<RANK> LDO;

  typename LDO::Ptr result = evolve(psi,sys,p);
  if (boost::shared_ptr<const SV> s = boost::dynamic_pointer_cast<const SV>(result)) {
    return arrayToNumpy<CArray<RANK>,RANK>(s->getArray());
  }
  else if(boost::shared_ptr<const DO> d = boost::dynamic_pointer_cast<const DO>(result)) {
    return arrayToNumpy<CArray<2*RANK>,2*RANK>(d->getArray());
  }
  else {
    PyErr_SetString(PyExc_NotImplementedError, "Result type of evolve not implemented in Python.");
    throw_error_already_set();
  }
  // this should be never reached
  return object();
}


void export_25_Evolution()
{
  import_array();
  numeric::array::set_module_and_type("numpy", "ndarray");
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
  def("evolve", &py_evolve<r>, "Wrapper of :core:`evolve` with RANK="#r);
BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_ADD(PYTHON_HALF_RANK,1), EVOLVE_INSTANTIATIONS, data)

}


} // pythonext