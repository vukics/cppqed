// Copyright András Vukics 2021–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Archive.h"
#include "Pars.h"
#include "Utility.h"

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
template <typename T, typename State, typename Time=double>
concept system = requires (T&& t, ConstReference<State> stateIn, Reference<State> stateOut, Time time) { t(stateIn,stateOut,time) ; };


template <typename State, typename Time=double>
using SystemFunctional = std::function<void(ConstReference<State>, Reference<State>, Time)>;

/// Embodies the concept defined at https://www.boost.org/doc/libs/1_80_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/concepts/controlled_stepper.html
/**
 * Lots of decays in concepts are ugly, however, cf.
 * [this suggestion](http://developercommunity.visualstudio.com/t/c-concepts-should-be-decayed/1268876)
 * [and this Q&A](http://stackoverflow.com/questions/74999063/when-if-ever-c-concepts-must-be-decayed-do-it-at-concept-definition-or-at)
 */
template <typename T, typename State, typename Time=double>
concept controlled_stepper =
  std::convertible_to<typename std::decay_t<T>::state_type,std::decay_t<State>> &&
  std::convertible_to<typename std::decay_t<T>::deriv_type,std::decay_t<State>> &&
  std::convertible_to<typename std::decay_t<T>::time_type,std::decay_t<Time>> &&
  requires (T&& t, Time deltaT, SystemFunctional<State,Time> sys, Time time, State&& stateInOut ) { t.try_step( sys , stateInOut , time , deltaT ); };


template <typename T, typename State, typename Time=double>
concept engine = adaptive_timestep_keeper<T> && intro_outro_streamer<T> &&
  requires ( T&& t, Time deltaT, std::ostream& logStream, SystemFunctional<State,Time> sys, Time& time, State&& stateInOut ) {
    { step(t,deltaT,logStream,sys,time,stateInOut) }; /*}  && 
  requires ( T&& t, Time deltaT, std::ostream& logStream, System sys, Time& time, ConstReference<State> stateIn, State&& stateOut ) {
    { step(t,deltaT,logStream,sys,time,stateIn,stateOut) };*/ };


using LogControl=std::bitset<2>;

/// Aggregate condensing parameters concerning adaptive ODE evolution in the style of a popl::OptionParser
/** This could be more finely grained by factoring out logControl, but the present solution should be sufficient for most puposes */
template <typename BASE=Empty>
struct Pars : BASE
{
  double
    epsRel, ///< relative precision of ODE stepping
    epsAbs; ///< absolute precision ”

  LogControl logControl{0};

  Pars(popl::OptionParser& op) : BASE{op}
  {
    addTitle(add(add(op,
      "epsAbs","ODE engine absolute precision",epsAbsDefault,&epsAbs),
      "epsRel","ODE engine relative precision",epsRelDefault,&epsRel),
      "ODE_Engine generic parameters");
  }
      
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
concept ode_logger = outro_streamer<T> && requires (T&& t, size_t failedStepsLast) { logDerivsCall(t); logStep(t,failedStepsLast); };

struct DefaultLogger { size_t nDerivsCalls=0, nSteps=0, nFailedSteps=0; };

void logDerivsCall(DefaultLogger& dl) {++dl.nDerivsCalls;}
void logStep(DefaultLogger& dl, size_t failedStepsLast) {++dl.nSteps; dl.nFailedSteps+=failedStepsLast;}

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
  
  Base(Time dtInit, LogControl lc, auto&&... args)
    : ces{MakeControlledErrorStepper<CES>::_(std::forward<decltype(args)>(args)...)},
      dtTry(dtInit), logControl(lc) {}

  template <typename State> requires controlled_stepper<CES,State,Time>
  auto tryStep(system<State,Time> auto sys, Time& time, State&& stateInOut) {return ces.stepper.try_step(sys,std::forward<State>(stateInOut),time,dtTry);}

  template <typename State> requires controlled_stepper<CES,State,Time>
  auto tryStep(system<State,Time> auto sys, Time& time, ConstReference<State> stateIn, State&& stateOut) {return ces.stepper.try_step(sys,stateIn,time,std::forward<State>(stateOut),dtTry);}
  
  ControlledErrorStepperWrapper<CES> ces;
  
  Time dtDid=0., dtTry;
  
  Logger logger;
  
  LogControl logControl;


  friend double getDtDid(const Base& b) {return b.dtDid;}

  friend std::ostream& streamIntro(const Base& b, std::ostream& os)
  {
    return b.ces.streamIntro(os<<"ODE engine: " << StepperDescriptor<CES> << ". ") << std::endl;
  }

  friend std::ostream& streamOutro(const Base& b, std::ostream& os) {return b.logControl[0] ? streamOutro(b.logger,os) : os;}

  /// The signature is such that it matches the signature of step in trajectories without the trailing parameters
  /**
  * The possibility of FSAL steppers (that can reuse the derivative calculated in the last stage of the previous step)
  * is not included, that would require more involved logic.
  */
  friend void step(Base& b, typename CES::time_type deltaT, std::ostream& logStream, auto sys, typename CES::time_type& time, auto&&... states)
  {
    using Time=typename CES::time_type;

    const Time timeBefore=time;
    Time nextDtTry=( fabs(deltaT)<fabs(b.dtTry/nextDtTryCorrectionFactor) ? b.dtTry : 0. );

    if ( ( deltaT<0 && b.dtTry>0 ) || ( deltaT>0 && b.dtTry<0 ) ) b.dtTry=-b.dtDid; // Stepping backward

    if (fabs(b.dtTry)>fabs(deltaT)) b.dtTry=deltaT;

    size_t nFailedSteps=0;

    for (;
          // wraps the sys functional in order that the number of calls to it can be logged
          b.tryStep( [&] (const auto& y, auto& dydt, double t) {
            logDerivsCall(b.logger);
            return sys(y,dydt,t);
          },time,std::forward<decltype(states)>(states)...)!=bno::success;
          ++nFailedSteps) ;

    b.dtDid=time-timeBefore;
    if (nextDtTry) b.dtTry=nextDtTry;

    logStep(b.logger,nFailedSteps);

    if (b.logControl[1]) logStream<<"Number of failed steps in this timestep: "<<nFailedSteps<<std::endl;
  }


  template <typename Archive>
  friend Archive& stateIO(Base& b, Archive& ar) {return SerializeControlledErrorStepper<CES>::_(ar,b.ces.stepper) & b.logger & b.dtDid & b.dtTry;}


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
    return ControlledErrorStepperWrapper<bno::controlled_runge_kutta<ErrorStepper>>{
      make_controlled(epsRel,epsAbs,ErrorStepper()),epsRel,epsAbs};
  }

  template <typename P>
  static auto _(const P& p) {return _(p.epsRel,p.epsAbs);}

  
};


} // ode


/// from this point on, Time is double everywhere
template <typename StateType>
using ODE_EngineBoost = ode::Base<boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54<StateType>>>;


namespace ode {

template <typename StateType, typename P>
ODE_EngineBoost<StateType> makeBoost(double dtInit, const P& p)
{
  return {dtInit,p.logControl,p.epsRel,p.epsAbs};
}



} // ode


} // cppqedutils

