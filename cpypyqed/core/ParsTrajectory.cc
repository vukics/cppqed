// Copyright Raimar Sandner 2012–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#include "Core.h"

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

PARS_GETTER_SETTER(double,  ParsRun, T)
PARS_GETTER_SETTER(int,     ParsRun, dc)
PARS_GETTER_SETTER(double,  ParsRun, Dt)
PARS_GETTER_SETTER(long,    ParsRun, NDt)
PARS_GETTER_SETTER(string,  ParsRun, ofn)
PARS_GETTER_SETTER(string,  ParsRun, initialFileName)
PARS_GETTER_SETTER(int,     ParsRun, precision)
PARS_GETTER_SETTER(bool,    ParsRun, displayInfo)
PARS_GETTER_SETTER(unsigned,ParsRun, sdf)


PARS_GETTER_SETTER(double, ParsEvolved, epsRel)
PARS_GETTER_SETTER(double, ParsEvolved, epsAbs)
PARS_GETTER_SETTER(evolved::SteppingFunction, ParsEvolved, sf)
PARS_GETTER_SETTER(double, ParsEvolved, nextDtTryCorrectionFactor)

void export_10_ParsTrajectory()
{
  {
    scope namespaceScope = evolvedNameSpace;
    enum_<evolved::SteppingFunction>("SF", "Wrapper of :core:`evolved::SteppingFunction`")
      .value("RK8PD",  evolved::SF_RK8PD)
      .value("RKCK",   evolved::SF_RKCK)
    ;
  }

  {
    scope namespaceScope = trajectoryNameSpace;
    class_<ParsRun>
      (
        "ParsRun",
        "Wrapper of :core:`trajectory::ParsRun`",
        init<parameters::ParameterTable&, optional<const string&> >()
          [with_custodian_and_ward<1,2>()]
      )
      .PARS_PROPERTY(T)
      .PARS_PROPERTY(dc)
      .PARS_PROPERTY(Dt)
      .PARS_PROPERTY(NDt)
      .PARS_PROPERTY(ofn)
      .PARS_PROPERTY(initialFileName)
      .PARS_PROPERTY(precision)
      .PARS_PROPERTY(displayInfo)
      .PARS_PROPERTY(sdf)
    ;
    class_<ParsEvolved>
      (
        "ParsEvolved",
        "Wrapper of :core:`trajectory::ParsEvolved`",
        init<parameters::ParameterTable&, optional<const string&> >()
          [with_custodian_and_ward<1,2>()]
      )
      .PARS_PROPERTY(epsRel)
      .PARS_PROPERTY(epsAbs)
      .PARS_PROPERTY(sf)
      .PARS_PROPERTY(nextDtTryCorrectionFactor)
    ;
  }
}

}