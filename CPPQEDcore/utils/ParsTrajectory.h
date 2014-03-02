/// \briefFile{Defines parameter aggergates to Trajectory.h}
// -*- C++ -*-
#ifndef CPPQEDCORE_UTILS_PARSTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_PARSTRAJECTORY_H_INCLUDED

#include "TrajectoryFwd.h"

#include "Evolved.h"
#include "FormDoubleFwd.h"
#include "Pars.h"


namespace trajectory {

/// Aggregate condensing parameters concerning adaptive ODE evolution (cf. Adaptive::Adaptive()) in the style of a parameters::ParameterTable
/**
 * The suit of ParsEvolved, ParsRun, ParsStochastic demonstrates a technique to facilitate reading parameters from command line and using them to initialize a (hierarchy of) class(es).
 * 
 * The same technique is to be used in elements when there is a detailed hierarchy of elements.
 * 
 * In principle the use of such `%Pars…` classes could be avoided by initializing classes with a parameters::ParameterTable. E.g.
 * 
 *     class X {
 * 
 *       double& y1;
 *       double& y2;
 * 
 *     public:
 *       X(Pars::ParameterTable& p) : y1(p.add("y1","y1 description",y1default)), y2(p.add("y2","y2 description",y2default)) {}
 * 
 *     };
 * 
 * Avoiding this technique is a design decision: its drawback would be that too many things (virtually everything) would depend on the \link parameters parameter-bundle\endlink,
 * and every parameter inside classes has to be declared as reference, which then makes it difficult to initialize the class in any other way than with a parameters::ParameterTable.
 *
 */
struct ParsEvolved
{
  static const double epsRelDefault; ///< The ultimate default of \link ParsEvolved::epsRel epsRel\endlink in the framework
  static const double epsAbsDefault; ///< ” for \link PrasEvolved::epsAbs epsAbs\endlink

  double
    &epsRel, ///< relative precision of ODE stepping (cf. evolved::TimeStepBookkeeper)
    &epsAbs; ///< absolute precision ”

  evolved::SteppingFunction& sf; ///< \link evolved::SteppingFunction stepping-function type\endlink
  double &nextDtTryCorrectionFactor; ///< cf. evolved::MakerGSL::MakerGSL()
  
  /// All `%Pars…` classes are constructed taking a parameters::ParameterTable, to register the parameters on
  /** 
   * This occurs via the parameters::ParameterTable::addMod member function returning a reference to the registered parameter wherewith the public attributes
   * (like ParsEvolved::epsRel) get initialized.
   */
  ParsEvolved(parameters::ParameterTable&, ///<[in/out] the table to register the new parameters on
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
  long &NDt; ///< number of deltaT intervals in \link run(Trajectory&, long, double, unsigned, const std::string&, const std::string&, int, bool, bool) deltaT-mode\endlink
  std::string &ofn, &initialFileName;

  formdouble::Zero &precision; ///< the overall precision of trajectory display \see FormDouble::overallPrecision

  bool &displayInfo, &firstStateDisplay;

  unsigned &sdf;
  
  ParsRun(parameters::ParameterTable&, const std::string& mod="");

  /// Corresponds to parameters::ParameterTable::getParsedCommandLine
  const std::string getParsedCommandLine() const {return *parsedCommandLine_;}

private:
  const parameters::ParameterTable::ParsedCommandLine parsedCommandLine_;

};


/// An example of a simple unification of `%Pars…` classes, here ParsRun & ParsEvolved
struct Pars : ParsRun, ParsEvolved
{
  Pars(parameters::ParameterTable& p, const std::string& mod="") : ParsRun(p,mod), ParsEvolved(p,mod) {}
};


} // trajectory

#endif // CPPQEDCORE_UTILS_PARSTRAJECTORY_H_INCLUDED
