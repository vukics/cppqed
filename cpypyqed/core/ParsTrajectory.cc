// -*- C++ -*-
#include "PythonExtension.h"

#include "EvolvedGSL.h"
#include "FormDouble.h"
#include "Pars.h"
#include "ParsTrajectory.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;

using trajectory::ParsEvolved;
using trajectory::ParsRun;
using std::string;

namespace pythonext{

PARS_GETTER_SETTER(double, ParsRun, T)
PARS_GETTER_SETTER(int,    ParsRun, dc)
PARS_GETTER_SETTER(double, ParsRun, Dt)
PARS_GETTER_SETTER(long,   ParsRun, NDt)
PARS_GETTER_SETTER(string, ParsRun, ofn)

PARS_GETTER_SETTER(double, ParsEvolved, epsRel)
PARS_GETTER_SETTER(double, ParsEvolved, epsAbs)




PARS_GETTER_SETTER(int,    Pars, precision)
PARS_GETTER_SETTER(bool,   Pars, displayInfo)
PARS_GETTER_SETTER(evolved::SteppingFunction, Pars, sf)
PARS_GETTER_SETTER(double, Pars, nextDtTryCorrectionFactor)

void export_ParsTrajectory()
{
  enum_<evolved::SteppingFunction>("SteppingFunction")
    .value("SF_RK8PD",  evolved::SF_RK8PD)
    .value("SF_RKCK",   evolved::SF_RKCK)
  ;
  class_<Pars>
    (
      "ParsTrajectory",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(T)
    .PARS_PROPERTY(epsRel)
    .PARS_PROPERTY(epsAbs)
    .PARS_PROPERTY(dc)
    .PARS_PROPERTY(Dt)
    .PARS_PROPERTY(NDt)
    .PARS_PROPERTY(ofn)
    .PARS_PROPERTY(autoStop)
    .PARS_PROPERTY(precision)
    .PARS_PROPERTY(displayInfo)
    .PARS_PROPERTY(sf)
    .PARS_PROPERTY(nextDtTryCorrectionFactor)
  ;
}

}