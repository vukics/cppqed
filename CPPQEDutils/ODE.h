// Copyright András Vukics 2021. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_UTILS_ODE_H_INCLUDED
#define CPPQEDCORE_UTILS_ODE_H_INCLUDED

#include "MathExtensions.h"
#include "Pars.h"

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif // BZ_HAVE_BOOST_SERIALIZATION

#include <boost/numeric/odeint.hpp>

#include <functional>
#include <iostream>
#include <optional>


namespace cppqedutils {

/*
template <typename E>
concept ordinary_differential_equation_engine = requires(E e) {
  { e.step(); } -> size_t;
};
*/

namespace ode_engine {

inline static const double epsRelDefault=1e-6; ///< The ultimate default of \link ode_engine::Pars::epsRel epsRel\endlink in the framework
inline static const double epsAbsDefault=1e-12; ///< ” for \link ode_engine::Pars::epsAbs epsAbs\endlink


template <typename Time, typename State>
using SystemFunctional_t=std::function<void(const State&, State&, Time)>;


/// Aggregate condensing parameters concerning adaptive ODE evolution in the style of a parameters::Table
/** If necessary, it can be made customizable by an additional template parameter, but a very sensible default can be provided */
template <typename BASE = parameters::Empty>
struct Pars : BASE
{
  double
    &epsRel, ///< relative precision of ODE stepping
    &epsAbs; ///< absolute precision ”

  int &logLevel; ///< governs how much logging information is streamed during a Trajectory run
  
  /// All `%Pars…` classes are constructed taking a parameters::Table, to register the parameters on
  /** 
   * This occurs via the parameters::Table::add member function returning a reference to the registered parameter wherewith the public attributes
   * (like Pars::epsRel) get initialized.
   */
  Pars(parameters::Table& p, ///<[in/out] the table to register the new parameters on
       const std::string& mod=""    ///<[in] possible modifier suffix
       )
    : BASE{p,mod},
      epsRel(p.addTitle("ODE_Engine generic parameters",mod).add("eps",mod,"ODE stepper relative precision",epsRelDefault)),
      epsAbs(p.add("epsAbs",mod,"ODE stepper absolute precision",epsAbsDefault)),
      logLevel(p.add("logLevel",mod,"logging level",0)) {}
      
};


namespace bno=boost::numeric::odeint;


template <typename ControlledErrorStepper>
struct ControlledErrorStepperParameters
{
  const double epsRel, epsAbs;
  
  std::ostream& stream(std::ostream& os) const {return os << "Parameters: epsRel=" << epsRel << " epsAbs=" << epsAbs;}
};


template <typename Stepper>
constexpr auto StepperDescriptor=std::nullopt;


template <typename ControlledErrorStepper>
struct MakeControlledErrorStepper;


template <typename ControlledErrorStepper>
struct SerializeControlledErrorStepper
{
  template <typename Archive>
  static Archive& _(Archive& ar, ControlledErrorStepper)
  {
    return ar;
  }
};


class DefaultLogger
{
public:
  void logDerivsCall() {++nDerivsCalls_;}
  void logStep() {++nSteps_;}
  void logFailedSteps(size_t failedStepsLast) {nFailedSteps_+=failedStepsLast;}
  
  std::ostream& logOnEnd(std::ostream& os) const
  {
    return os<<"\nTotal number of ODE steps: "<<nSteps_
             <<"\nNumber of failed ODE steps: "<<nFailedSteps_
             <<"\nNumber of calls to function calculating RHS for ODE: "<<nDerivsCalls_<<std::endl;
  }

private:
#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & nDerivsCalls_ & nSteps_ & nFailedSteps_;}
#endif // BZ_HAVE_BOOST_SERIALIZATION

  size_t nDerivsCalls_=0, nSteps_=0, nFailedSteps_=0;

};


template <typename ControlledErrorStepper, typename Logger=DefaultLogger>
class Base
{
public:
  /**
   * \param nextDtTryCorrectionFactor
   * Since dtTry is at most deltaT, dtTry will drop severely if deltaT is very small (e.g. when completing a given interval at its end).
   * Since this drop can be several orders of magnitude, this must not be allowed as it would slow down the simulation extremely.
   * To patch this, we check before the actual timestep whether we are in the case of a *very small* `deltaT` (`< dtTry/nextDtTryCorrectionFactor`).
   * In such a case the actual dtTry is cached and used unchanged in the next step. That is, such a tiny step is excluded from stepsize control.
   */
  inline static const double nextDtTryCorrectionFactor=10.;
  
  
  using Time=typename ControlledErrorStepper::time_type;
  
  template <typename ... Args>
  Base(Time dtInit, int logLevel, Args&&... args)
    : ces_{MakeControlledErrorStepper<ControlledErrorStepper>::_(std::forward<Args>(args)...)},
      dtTry_(dtInit), logLevel_(logLevel) {}

  template <typename ParsBase>
  Base(Time dtInit, const Pars<ParsBase>& p)
    : ces_{MakeControlledErrorStepper<ControlledErrorStepper>::_(p)},
      dtTry_(dtInit), logLevel_(p.logLevel) {}

  /// The signature is such that it matches the signature of step in trajectories without the trailing parameters
  template <typename System, typename ... States>
  void step(Time deltaT, std::ostream& logStream, System sys, Time& time, States&&... states);

  Time getDtDid() const {return dtDid_;}

  Time getDtTry() const {return dtTry_;}
  
  void setDtTry(Time dtTry) {dtTry_=dtTry;}

  std::ostream& streamParameters(std::ostream& os) const {
    return std::get<1>(ces_).stream(os<<"ODE engine implementation: " << StepperDescriptor<ControlledErrorStepper> << ". ")<< std::endl;
  }

  template <typename Archive>
  Archive& stateIO(Archive& ar) {return SerializeControlledErrorStepper<ControlledErrorStepper>::_(ar,std::get<0>(ces_)) & logger_ & dtDid_ & dtTry_;}
  
  std::ostream& logOnEnd(std::ostream& os) const {return logLevel_ > 0 ? logger_.logOnEnd(os) : os;}
  
private:
  template <typename System, typename State>
  auto tryStep(System sys, Time& time, State&& stateInOut) {return std::get<0>(ces_).try_step(sys,std::forward<State>(stateInOut),time,dtTry_);}

  template <typename System, typename State>
  auto tryStep(System sys, Time& time, const State& stateIn, State&& stateOut) {return std::get<0>(ces_).try_step(sys,stateIn,time,std::forward<State>(stateOut),dtTry_);}
  
  std::tuple<ControlledErrorStepper,
             ControlledErrorStepperParameters<ControlledErrorStepper>> ces_;
  
  Time dtDid_=0., dtTry_;
  
  Logger logger_;
  
  int logLevel_;

};


/// Specializations for Boost.Odeint controlled_runge_kutta & runge_kutta_cash_karp54

template <typename ErrorStepper>
inline const std::string StepperDescriptor<bno::controlled_runge_kutta<ErrorStepper>> = "Boost.Odeint controlled stepper " + StepperDescriptor<ErrorStepper>;


template <typename StateType>
inline const std::string StepperDescriptor<bno::runge_kutta_cash_karp54<StateType>> = "RKCK54";


template <typename ErrorStepper>
struct MakeControlledErrorStepper<bno::controlled_runge_kutta<ErrorStepper>>
{
  static auto _(double epsRel, double epsAbs)
  {
    return std::make_tuple(make_controlled(epsRel,epsAbs,ErrorStepper()),
                           ControlledErrorStepperParameters<bno::controlled_runge_kutta<ErrorStepper>>{epsRel,epsAbs});
  }

  template <typename ParsBase>
  static auto _(const Pars<ParsBase>& p)
  {
    return _(p.epsRel,p.epsAbs);
  }

  
};


} // ode_engine


/// from this point on, every Time is double
template <typename StateType>
using ODE_EngineBoost = ode_engine::Base<boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54<StateType>>>;



} // cppqedutils


/**
 * The possibility of FSAL steppers (that can reuse the derivative calculated in the last stage of the previous step)
 * is not included, that would require more involved logic.
 */
template <typename ControlledErrorStepper, typename Logger> template <typename System, typename ... States>
void cppqedutils::ode_engine::Base<ControlledErrorStepper,Logger>::step(Time deltaT, std::ostream& logStream, System sys, Time& time, States&&... states)
{  
  const Time timeBefore=time;
  
  Time nextDtTry=( fabs(deltaT)<fabs(dtTry_/nextDtTryCorrectionFactor) ? dtTry_ : 0. );
  
  if (sign(deltaT)!=sign(dtTry_)) dtTry_=-dtDid_; // Stepping backward
  
  if (fabs(dtTry_)>fabs(deltaT)) dtTry_=deltaT;
  
  size_t nFailedSteps=0;
  
  for (;
        // wraps the sys functional in order that the number of calls to it can be logged
        tryStep([this,sys](const auto& y, auto& dydt, double t) {
          logger_.logDerivsCall();
          return sys(y,dydt,t);
        },time,std::forward<States>(states)...)!=bno::success;
        ++nFailedSteps) ;

  dtDid_=time-timeBefore;
  if (nextDtTry) dtTry_=nextDtTry;
  
  logger_.logStep();
  logger_.logFailedSteps(nFailedSteps);
  
  if (logLevel_>3) logStream<<"Number of failed steps in this timestep: "<<nFailedSteps<<std::endl;
}


#endif // CPPQEDCORE_UTILS_ODE_H_INCLUDED
