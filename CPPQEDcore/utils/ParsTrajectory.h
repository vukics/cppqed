// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines parameter aggergates to Trajectory.h}
#ifndef CPPQEDCORE_UTILS_PARSTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_PARSTRAJECTORY_H_INCLUDED

#include "Evolved.h"
#include "FormDouble.h"
#include "ParsFwd.h"


namespace trajectory {

/// Aggregate condensing parameters concerning adaptive ODE evolution (cf. Adaptive::Adaptive()) in the style of a parameters::Table
/**
 * The suit of ParsEvolved, ParsRun, ParsStochastic demonstrates a technique to facilitate reading parameters from command line and using them to initialize a (hierarchy of) class(es).
 * 
 * The same technique is to be used in elements when there is a detailed hierarchy of elements.
 * 
 * In principle the use of such `%Pars…` classes could be avoided by initializing classes with a parameters::Table. E.g.
 * 
 *     class X {
 * 
 *       double& y1;
 *       double& y2;
 * 
 *     public:
 *       X(parameters::Table& p) : y1(p.add("y1","y1 description",y1default)), y2(p.add("y2","y2 description",y2default)) {}
 * 
 *     };
 * 
 * Avoiding this technique is a design decision: its drawback would be that too many things (virtually everything) would depend on the \link parameters parameter-bundle\endlink,
 * and every parameter inside classes has to be declared as reference, which then makes it difficult to initialize the class in any other way than with a parameters::Table.
 *
 */
struct ParsEvolved
{
  static const double epsRelDefault; ///< The ultimate default of \link ParsEvolved::epsRel epsRel\endlink in the framework
  static const double epsAbsDefault; ///< ” for \link ParsEvolved::epsAbs epsAbs\endlink

  double
    &epsRel, ///< relative precision of ODE stepping (cf. evolved::TimeStepBookkeeper)
    &epsAbs; ///< absolute precision ”

  evolved::SteppingFunction& sf; ///< \link evolved::SteppingFunction stepping-function type\endlink
  double &nextDtTryCorrectionFactor; ///< cf. evolved::MakerGSL::MakerGSL()
  
  int &logLevel; ///< governs how much logging information is displayed during a Trajectory run \see quantumtrajectory::mcwf::Logger for a usecase
  
  /// All `%Pars…` classes are constructed taking a parameters::Table, to register the parameters on
  /** 
   * This occurs via the parameters::Table::addMod member function returning a reference to the registered parameter wherewith the public attributes
   * (like ParsEvolved::epsRel) get initialized.
   */
  ParsEvolved(parameters::Table&, ///<[in/out] the table to register the new parameters on
              const std::string& mod=""    ///<[in] possible modifier suffix
             );
  
};


/// Parameters corresponding to the different versions of run()
/** \see ParsEvolved for general explanation of the workings of `%Pars…` classes */
struct ParsRun
{
  double &T; ///< endtime of the run
  int &dc;
  double &Dt;
  long &NDt; ///< number of deltaT intervals in \link trajectory::Trajectory::run deltaT-mode\endlink
  std::string &ofn, &initialFileName;

  formdouble::Zero &precision; ///< the overall precision of trajectory display \see FormDouble::overallPrecision

  bool &displayInfo, &firstStateDisplay;

  unsigned &sdf;

  double &autoStopEpsilon;

  unsigned &autoStopRepetition;
  
  ParsRun(parameters::Table&, const std::string& mod="");

  /// Corresponds to parameters::Table::getParsedCommandLine
  const std::string getParsedCommandLine() const {return *parsedCommandLine_;}

private:
  const parameters::Table::ParsedCommandLine parsedCommandLine_;

};


/// An example of a simple unification of `%Pars…` classes, here ParsRun & ParsEvolved
struct Pars : ParsRun, ParsEvolved
{
  Pars(parameters::Table& p, const std::string& mod="") : ParsRun(p,mod), ParsEvolved(p,mod) {}
};


} // trajectory

#endif // CPPQEDCORE_UTILS_PARSTRAJECTORY_H_INCLUDED
