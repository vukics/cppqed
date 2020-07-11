// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "Core.h"
#include "blitz2numpy.h"

#include "ParsPropertyMacros.h"

#include "EvolvedGSL.tcc"
#include "Pars.h"

#include "QuantumSystem.h"

#include "DensityOperator.h"
#include "StateVector.h"
#include "LazyDensityOperator.h"

#include "EnsembleMCWF.tcc"
#include "Evolution.tcc"
#include "Master.tcc"
#include "ParsMCWF_Trajectory.h"

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <boost/python/import.hpp>

using namespace boost::python;
using namespace evolution;

using trajectory::ParsRun;
using quantumdata::StateVector;
using structure::QuantumSystem;

namespace pythonext {

PARS_GETTER_SETTER(Method, Pars<>, evol)
PARS_GETTER_SETTER(bool, Pars<>, negativity)

template<int RANK>
object py_evolve(const numpy::ndarray &array,
                 typename structure::QuantumSystem<RANK>::Ptr sys,
                 const evolution::Pars<>& p)
{
  using namespace quantumdata;

  object qs = import("cpypyqed.tools.quantumstate");
  object py_SV = qs.attr("StateVector");
  object py_DO = qs.attr("DensityOperator");

  typedef StateVector        <RANK> SV;
  typedef DensityOperator    <RANK> DO;
  typedef LazyDensityOperator<RANK> LDO;

  typename LDO::Ptr result;

  if (array.get_nd() == RANK) {
    StateVector<RANK> psi(numpyToArray<dcomp,RANK>(array),quantumdata::ByReference());
    result = evolve(psi,sys,p);
  }
  else if (array.get_nd() == 2*RANK) {
    DensityOperator<RANK> rho(numpyToArray<dcomp,2*RANK>(array),quantumdata::ByReference());
    result = evolve(rho,sys,p);
  }
  else {
    PyErr_SetString(PyExc_ValueError, (std::string("Expected array of rank ") + std::to_string(RANK) + std::string(" or ") + std::to_string(2*RANK) + std::string(".")).c_str());
    throw_error_already_set();
  }

  if (std::shared_ptr<const SV> s = std::dynamic_pointer_cast<const SV>(result)) {
    return py_SV(*s);
  }
  else if(std::shared_ptr<const DO> d = std::dynamic_pointer_cast<const DO>(result)) {
    return py_DO(*d);
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
  {
    scope namespaceScope = evolutionNameSpace;
    enum_<Method>("Method", "Wrapper of :core:`evolution::Method`")
      .value("SINGLE",   SINGLE)
      .value("ENSEMBLE", ENSEMBLE)
      .value("MASTER",   MASTER)
    ;
    class_<Pars<>, bases<ParsRun,quantumtrajectory::mcwf::Pars> >
      (
        "Pars", "Wrapper of :core:`evolution::Pars`",
      init<parameters::Table&, optional<const std::string&> >()
      [with_custodian_and_ward<1,2>()]
      )
      .PARS_PROPERTY(evol)
      .PARS_PROPERTY(negativity)
    ;
  }


  // In the docstring of these wrapper functions: working around a bug in doxylink which strips 'struct' from the 'structure' namespace
#define EVOLVE_INSTANTIATIONS(z,r,data) \
  def("evolve", (object(*)(const numpy::ndarray &array,typename QuantumSystem<r>::Ptr,const evolution::Pars<>&))&py_evolve<r>, "Wrapper of :core:`evolve <Evolution_.h::evolve(quantumdata::StateVector< RANK >&, const ure::QuantumSystem< RANK >&, const evolution::Pars<>&)>` with RANK="#r ".");
BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_ADD(PYTHON_HALF_RANK,1), EVOLVE_INSTANTIATIONS, data)

}


} // pythonext
