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


/*
template <typename E>
concept ordinary_differential_equation_engine = requires(E e) {
  { e.step(); }
};
*/

namespace ode_engine {

namespace bno=boost::numeric::odeint;

template <typename ControlledErrorStepper>
struct ControlledErrorStepperParameters;

template <typename ControlledErrorStepper>
constexpr auto ControlledErrorStepperDescriptor=std::nullopt;

template <typename ErrorStepper>
constexpr auto ErrorStepperDescriptor=std::nullopt;

template <typename ErrorStepper>
auto makeControlledErrorStepper(bno::controlled_runge_kutta<ErrorStepper>, double epsRel, double epsAbs)
{
  return std::make_tuple(make_controlled(epsRel,epsAbs,ErrorStepper()),
                         ControlledErrorStepperParameters<bno::controlled_runge_kutta<ErrorStepper>>{epsRel,epsAbs});
}


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
class ODE_EngineBase
{
public:
  using Time=typename ControlledErrorStepper::time_type;
  
  template <typename ... Args>
  ODE_EngineBase(Time dtInit, Args&&... args)
    : ces_{makeControlledErrorStepper(ControlledErrorStepper(),std::forward<Args>(args)...)},
      dtTry_(dtInit) {}
  
  template <typename System, typename ... States>
  size_t step(System sys, Time& time, Time deltaT, States&&... states) {
    if (mathutils::sign(deltaT)!=mathutils::sign(EvolvedIO<A>::getDtTry())) setDtTry(-EvolvedIO<A>::getDtDid()); // Stepping backward
    
    Time timeBefore=time;
    dtTry_=std::min(dtTry_,deltaT);
    
    size_t nFailedSteps=0;
    
    for (;
         // wraps the sys functional in order that the number of calls to it can be logged
         stepImpl([this,sys](const auto& y, auto& dydt, double t) {
           logger_.logDerivsCall();
           return sys(y,dydt,t);
         },time,std::forward<States>(states)...)!=bno::success;
         ++nFailedSteps) ;

    dtDid_=time-timeBefore;
    
    logger_.logStep();
    logger_.logFailedSteps(nFailedSteps);
    
    return nFailedSteps;
  }

  Time getDtDid() const {return dtDid_;}

  Time getDtTry() const {return dtTry_;}
  
  void setDtTry(Time dtTry) {dtTry_=dtTry;}

  std::ostream& streamParameters(std::ostream& os) const {
    return os<<"ODE engine implementation: " << ControlledErrorStepperDescriptor<ControlledErrorStepper> << " " << std::get<1>(ces_) << std::endl;
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

  std::tuple<ControlledErrorStepper,ControlledErrorStepperParameters<ControlledErrorStepper>> ces_;
  
  Time dtDid_=0., dtTry_;
  
  Logger logger_;

};


template <typename ErrorStepper>
struct ControlledErrorStepperParameters<bno::controlled_runge_kutta<ErrorStepper>>
{
  const double epsRel, epsAbs;
};


template <typename ErrorStepper>
constexpr auto ControlledErrorStepperDescriptor<bno::controlled_runge_kutta<ErrorStepper>> = "Boost.Odeint controlled stepper " + ErrorStepperDescriptor<ErrorStepper>;

template <typename StateType>
constexpr auto ErrorStepperDescriptor<bno::runge_kutta_cash_karp54<StateType>> = "RKCK54";

} // ode_engine


template <typename StateType>
using ODE_EngineBoost = ode_engine::ODE_EngineBase<boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54<StateType>>>;


namespace cpputils {


/// \name Generic evolution functions
//@{
/// evolves for exactly time `deltaT`
/** \tparam E type of the object to evolve. Implicit interface assumed: member function named step with signature `...(double)` */
template<typename E>
void evolve(E& e, double deltaT, std::ostream& logStream=std::clog)
{
  double endTime=e.getTime()+deltaT;
  while (double dt=endTime-e.getTime()) e.step(dt,logStream);
}


/// evolves up to exactly time `t` \copydetails evolve
template<typename E>
void evolveTo(E& e, double t, std::ostream& logStream=std::clog)
{
  evolve(e,t-e.getTime(),logStream);
}
//@}


} // cpputils


#endif // CPPQEDCORE_UTILS_ODE_H_INCLUDED
