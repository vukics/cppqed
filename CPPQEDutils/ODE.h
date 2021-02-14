// Copyright András Vukics 2021. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_UTILS_ODE_H_INCLUDED
#define CPPQEDCORE_UTILS_ODE_H_INCLUDED

#include "MathExtensions.h"

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif // BZ_HAVE_BOOST_SERIALIZATION

#include <boost/numeric/odeint.hpp>

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


namespace bno=boost::numeric::odeint;

template <typename ControlledErrorStepper>
struct ControlledErrorStepperParameters;

template <typename Stepper>
constexpr auto StepperDescriptor=std::nullopt;


template <typename ControlledErrorStepper>
struct MakeControlledErrorStepper;


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
             <<"\nNumber of calls of function calculating RHS for ODE: "<<nDerivsCalls_<<std::endl;
  }

private:
#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & nDerivsCalls_ & nSteps_ & nFailedSteps_;}
#endif // BZ_HAVE_BOOST_SERIALIZATION

  size_t nDerivsCalls_, nSteps_, nFailedSteps_;

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
   * In such a case the present dtTry is cached and used unchanged in the next step. That is, such a tiny step is excluded from stepsize control.
   */
  inline static const double nextDtTryCorrectionFactor=50.;
  
  
  using Time=typename ControlledErrorStepper::time_type;
  
  template <typename ... Args>
  Base(Time dtInit, Args&&... args)
    : ces_{MakeControlledErrorStepper<ControlledErrorStepper>{}(std::forward<Args>(args)...)},
      dtTry_(dtInit) {}
  
  /// The signature is such that it matches the signature of step in trajectories without the trailing parameters
  template <typename System, typename ... States>
  size_t step(Time deltaT, std::ostream& logStream, System sys, Time& time, States&&... states) {
    
    Time 
      timeBefore=time,
      nextDtTry=( fabs(deltaT)<fabs(dtTry_/nextDtTryCorrectionFactor) ? dtTry_ : 0. );
    
    if (sign(deltaT)!=sign(dtTry_)) dtTry_=-dtDid_; // Stepping backward
    
    if (fabs(dtTry_)>fabs(deltaT)) dtTry_=deltaT;
    
    size_t nFailedSteps=0;
    
    for (;
         // wraps the sys functional in order that the number of calls to it can be logged
         stepImpl([this,sys](const auto& y, auto& dydt, double t) {
           logger_.logDerivsCall();
           return sys(y,dydt,t);
         },time,std::forward<States>(states)...)!=bno::success;
         ++nFailedSteps) ;

    dtDid_=time-timeBefore;
    if (nextDtTry) dtTry_=nextDtTry;
    
    logger_.logStep();
    logger_.logFailedSteps(nFailedSteps);
    
    return nFailedSteps;
  }

  Time getDtDid() const {return dtDid_;}

  Time getDtTry() const {return dtTry_;}
  
  void setDtTry(Time dtTry) {dtTry_=dtTry;}

  std::ostream& streamParameters(std::ostream& os) const {
    return std::get<1>(ces_).stream(os<<"ODE engine implementation: " << StepperDescriptor<ControlledErrorStepper> << ". ")<< std::endl;
  }
  
private:
  template <typename System, typename State>
  auto stepImpl(System sys, Time& time, State&& stateInOut) {return std::get<0>(ces_).try_step(sys,std::forward<State>(stateInOut),time,dtTry_);}

  template <typename System, typename State>
  auto stepImpl(System sys, Time& time, const State& stateIn, State&& stateOut) {return std::get<0>(ces_).try_step(sys,stateIn,time,std::forward<State>(stateOut),dtTry_);}
  
#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & std::get<0>(ces_) & logger_ & dtDid_ & dtTry_;}
#endif // BZ_HAVE_BOOST_SERIALIZATION

  std::tuple<ControlledErrorStepper,
             ControlledErrorStepperParameters<ControlledErrorStepper>> ces_;
  
  Time dtDid_=0., dtTry_;
  
  Logger logger_;

};


/// Specializations for Boost.Odeint controlled_runge_kutta & runge_kutta_cash_karp54

template <typename ErrorStepper>
struct ControlledErrorStepperParameters<bno::controlled_runge_kutta<ErrorStepper>>
{
  const double epsRel, epsAbs;
  
  std::ostream& stream(std::ostream& os) const {return os << "Parameters: epsRel=" << epsRel << " epsAbs=" << epsAbs;}
  
};


template <typename ErrorStepper>
inline const std::string StepperDescriptor<bno::controlled_runge_kutta<ErrorStepper>> = "Boost.Odeint controlled stepper " + StepperDescriptor<ErrorStepper>;


template <typename StateType>
inline const std::string StepperDescriptor<bno::runge_kutta_cash_karp54<StateType>> = "RKCK54";


template <typename ErrorStepper>
struct MakeControlledErrorStepper<bno::controlled_runge_kutta<ErrorStepper>>
{
  auto operator()(double epsRel, double epsAbs) {
    return std::make_tuple(make_controlled(epsRel,epsAbs,ErrorStepper()),
                           ControlledErrorStepperParameters<bno::controlled_runge_kutta<ErrorStepper>>{epsRel,epsAbs});
  }
};


} // ode_engine


/// from this point on, every Time is double
template <typename StateType>
using ODE_EngineBoost = ode_engine::Base<boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54<StateType>>>;



} // cppqedutils


#endif // CPPQEDCORE_UTILS_ODE_H_INCLUDED
