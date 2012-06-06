// -*- C++ -*-
#ifndef _PARS_TRAJECTORY_H
#define _PARS_TRAJECTORY_H

#include "TrajectoryFwd.h"

#include "EvolvedGSL.h"
#include "ParsFwd.h"


namespace trajectory {

/*
  The following suit of Pars(Stochastic)Trajectory, + eg ParsMCWF_Trajectory in C++QED demonstrates a technique to facilitate reading parameters from command line and using them to initialize a (hierarchy of) class(es).

  The same technique is to be used in quantum/elements when there is a detailed hierarchy of elements.

  In principle the use of Pars(...) classes could be avoided by initializing classes with a ParameterTable. Eg

  class X {
  double& y1;
  double& y2;

  public:
  X(Pars::ParameterTable& p) : 
  y1(p.add("y1","y1 description",y1default)), y2(p.add("y2","y2 description",y2default)) {}

  };

  Avoiding this technique is a design decision: the drawback is that too many things (virtually everything) would depend on Pars, and every parameter inside classes has to be declared as reference, which then makes it difficult to initialize the class NOT with a ParameterTable.

*/

struct ParsTrajectory {

  static const double epsRelDefault;
  static const double epsAbsDefault;

  double &T, &epsRel, &epsAbs;
  int &dc;
  double &Dt;
  std::string &ofn;
  double &autoStop;
  int &precision;
  bool &displayInfo;

  evolved::SteppingFunction& sf;
  double &nextDtTryCorretionFactor;

  ParsTrajectory(parameters::ParameterTable&, const std::string& mod="");

  virtual ~ParsTrajectory() {}

};


} // trajectory

#endif // _PARS_TRAJECTORY_H
