// Copyright András Vukics 2021–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Archive.h"
#include "Pars.h"

#include <boost/numeric/odeint.hpp>

#include <bitset>
#include <iostream>
#include <optional>


namespace cppqedutils {


template <typename T> concept adaptive_timestep_keeper = requires (const T& t) { 
  { getDtDid(t) } -> std::convertible_to<double>; };

template <typename T> concept intro_streamer = requires ( const T& t, std::ostream& os ) { 
  { streamIntro(t,os) } -> std::convertible_to<std::ostream&>; };

template <typename T> concept outro_streamer = requires ( const T& t, std::ostream& os ) { 
  { streamOutro(t,os) } -> std::convertible_to<std::ostream&>; };

template <typename T> concept intro_outro_streamer = intro_streamer<T> && outro_streamer<T>;


namespace ode {

inline static const double epsRelDefault=1e-6; ///< The ultimate default of \link ode_engine::Pars::epsRel epsRel\endlink in the framework
inline static const double epsAbsDefault=1e-12; ///< ” for \link ode_engine::Pars::epsAbs epsAbs\endlink

/**
 * Since dtTry is at most deltaT, dtTry will drop severely if deltaT is very small (e.g. when completing a given interval at its end).
 * Since this drop can be several orders of magnitude, this must not be allowed as it would slow down the simulation extremely.
 * To patch this, we check before the actual timestep whether we are in the case of a *very small* `deltaT` (`< dtTry/nextDtTryCorrectionFactor`).
 * In such a case the actual dtTry is cached and used unchanged in the next step. That is, such a tiny step is excluded from stepsize control.
 */
inline static const double nextDtTryCorrectionFactor=10.;


/// Embodies the concept defined at https://www.boost.org/doc/libs/1_80_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/concepts/system.html
template <typename T, typename Time, typename State>
concept system = requires (T&& t, const State& stateIn, State& stateOut, Time time) { t(stateIn,stateOut,time) ; };

/// Embodies the concept defined at https://www.boost.org/doc/libs/1_80_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/concepts/controlled_stepper.html
template <typename T, typename Time, typename State, typename System>
concept controlled_stepper =
  system<System,Time,State> &&
  std::convertible_to<typename T::state_type,State> && std::convertible_to<typename T::deriv_type,System> && std::convertible_to<typename T::time_type,Time> &&
  requires (T&& t, System sys, Time time, Time dt, State&& x ) { t.try_step( sys , x , time , dt ); };


template <typename T, typename Time, typename State, typename System>
concept engine = 
  system<System,Time,State> && adaptive_timestep_keeper<T> && intro_outro_streamer<T> &&
  requires ( T&& t, Time deltaT, std::ostream& logStream, System sys, Time& time, State&& stateInOut ) { 
    { step(t,deltaT,logStream,sys,time,stateInOut) }; /*}  && 
  requires ( T&& t, Time deltaT, std::ostream& logStream, System sys, Time& time, const State& stateIn, State&& stateOut ) { 
    { step(t,deltaT,logStream,sys,time,stateIn,stateOut) };*/ };

 
/// Aggregate condensing parameters concerning adaptive ODE evolution in the style of a parameters::Table
/** If necessary, it can be made customizable by an additional template parameter, but a very sensible default can be provided */
struct Pars
{
  double
    epsRel, ///< relative precision of ODE stepping
    epsAbs; ///< absolute precision ”

  Pars(popl::OptionParser& op)
  {
    op.add<popl::Value<double>>("","epsRel","ODE engine relative precision",epsRelDefault,&epsRel);
    op.add<popl::Value<double>>("","epsAbs","ODE engine absolute precision",epsAbsDefault,&epsAbs);
  }
   /* : epsRel(p.addTitle("ODE_Engine generic parameters",mod).add("eps",mod,"ODE stepper relative precision",epsRelDefault)),
      epsAbs(p.add("epsAbs",mod,"ODE stepper absolute precision",epsAbsDefault)),
      lc(p.add("lc",mod,"logging level",0)) {}*/
      
};


namespace bno=boost::numeric::odeint;


/// Binds together a ControlledErrorStepper with its construction parameters, that unfortunalety cannot be recovered from the Boost steppers otherwise.
template <typename ControlledErrorStepper>
struct ControlledErrorStepperWrapper
{
  std::ostream& streamIntro(std::ostream& os) const {return os << "Parameters: epsRel=" << epsRel << ", epsAbs=" << epsAbs;}

  ControlledErrorStepper stepper;
  double epsRel, epsAbs;
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


template <typename T>
concept ode_logger = outro_streamer<T> && requires (T t) { logDerivsCall(t); logStep(t); logFailedSteps(t); };

struct DefaultLogger { size_t nDerivsCalls=0, nSteps=0, nFailedSteps=0; };

void logDerivsCall(DefaultLogger& dl) {++dl.nDerivsCalls;}
void logStep(DefaultLogger& dl) {++dl.nSteps;}
void logFailedSteps(DefaultLogger& dl, size_t failedStepsLast) {dl.nFailedSteps+=failedStepsLast;}

template<typename Archive>
inline void serialize(Archive & ar, DefaultLogger& dl, const unsigned int /*file_version*/) { ar & dl.nDerivsCalls & dl.nSteps & dl.nFailedSteps; }

std::ostream& streamOutro(const DefaultLogger& dl, std::ostream& os)
{
  return os<<"\nODE DefaultLogger\nTotal number of ODE steps: "<<dl.nSteps
           <<"\nNumber of failed ODE steps: "<<dl.nFailedSteps
           <<"\nNumber of calls to function calculating RHS for ODE: "<<dl.nDerivsCalls<<std::endl;
}


template <typename CES, 
//< here we cannot use the above concept, since at this point we don’t know the time, state, and system types – this comes only in the step function
          ode_logger Logger=DefaultLogger>
struct Base
{
  using Time=typename CES::time_type;
  
  using LogControl=std::bitset<2>;
  
  Base(Time dtInit, LogControl lc, auto&&... args)
    : ces{MakeControlledErrorStepper<CES>::_(std::forward<decltype(args)>(args)...)},
      dtTry(dtInit), logControl(lc) {}

  template <typename State, system<Time,State> SYS> requires controlled_stepper<CES,Time,State,SYS>
  auto tryStep(SYS sys, Time& time, State&& stateInOut) {return ces.stepper.try_step(sys,std::forward<State>(stateInOut),time,dtTry);}

  template <typename State, system<Time,State> SYS> requires controlled_stepper<CES,Time,State,SYS>
  auto tryStep(SYS sys, Time& time, const State& stateIn, State&& stateOut) {return ces.stepper.try_step(sys,stateIn,time,std::forward<State>(stateOut),dtTry);}
  
  ControlledErrorStepperWrapper<CES> ces;
  
  Time dtDid=0., dtTry;
  
  Logger logger;
  
  LogControl logControl;

  template <typename Archive>
  Archive& stateIO(Archive& ar) {return SerializeControlledErrorStepper<CES>::_(ar,ces.stepper) & logger & dtDid & dtTry;}

};


template <typename CES, typename Logger> auto getDtDid(const Base<CES,Logger>& b) {return b.dtDid;}

template <typename CES, typename Logger> std::ostream& streamIntro(const Base<CES,Logger>& b, std::ostream& os)
{
  return b.ces.streamIntro(os<<"ODE engine: " << StepperDescriptor<CES> << ". ") << std::endl;
}

template <typename CES, typename Logger> std::ostream& streamOutro(const Base<CES,Logger>& b, std::ostream& os)
{
  return b.logControl[0] ? streamOutro(b.logger,os) : os;
}


/// The signature is such that it matches the signature of step in trajectories without the trailing parameters
/**
 * The possibility of FSAL steppers (that can reuse the derivative calculated in the last stage of the previous step)
 * is not included, that would require more involved logic.
 */
template <typename CES, typename Logger>
void step(Base<CES,Logger>& b, typename CES::time_type deltaT, std::ostream& logStream, auto sys, typename CES::time_type& time, auto&&... states)
{
  using Time=typename CES::time_type;
  
  const Time timeBefore=time;
  Time nextDtTry=( fabs(deltaT)<fabs(b.dtTry/b.nextDtTryCorrectionFactor) ? b.dtTry : 0. );
  
  if (sign(deltaT)!=sign(b.dtTry)) b.dtTry=-b.dtDid; // Stepping backward
  
  if (fabs(b.dtTry)>fabs(deltaT)) b.dtTry=deltaT;
  
  size_t nFailedSteps=0;
  
  for (;
        // wraps the sys functional in order that the number of calls to it can be logged
        tryStep( [&] (const auto& y, auto& dydt, double t) {
          b.logger.logDerivsCall();
          return sys(y,dydt,t);
        },time,std::forward<decltype(states)>(states)...)!=bno::success;
        ++nFailedSteps) ;

  b.dtDid=time-timeBefore;
  if (nextDtTry) b.dtTry=nextDtTry;
  
  b.logger.logStep();
  b.logger.logFailedSteps(nFailedSteps);
  
  if (b.logControl[1]) logStream<<"Number of failed steps in this timestep: "<<nFailedSteps<<std::endl;
}


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
    return ControlledErrorStepperWrapper<bno::controlled_runge_kutta<ErrorStepper>>{
      make_controlled(epsRel,epsAbs,ErrorStepper()),epsRel,epsAbs};
  }

  static auto _(const Pars& p) {return _(p.epsRel,p.epsAbs);}

  
};


} // ode


/// from this point on, every Time is double
template <typename StateType>
using ODE_EngineBoost = ode::Base<boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54<StateType>>>;



} // cppqedutils

