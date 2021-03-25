// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED
#define CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED

#include "ArrayTraits.h"
#include "KeyPrinter.h"
#include "Trajectory.h"

#include <typeinfo>

namespace cppqedutils {

/// Class fully implementing the Adaptive interface by streaming (and serializing) the whole content of the evolved array
/**
 * Meant for all cases when simple ODE evolution is needed with intermittent streamings
 * 
 * <b>Example usage:</b> simulation of a complex driven damped harmonic oscillator mode described by the ODE \f[\ddot{y}+2\gamma\,\dot{y}+y=e^{i\,\omega t},\f]
 * where \f$\gamma\f$ is the damping rate and \f$\omega\f$ the driving frequency, and the timescale has been chosen such that the eigenfrequency is 1.
 * 
 * \include HarmonicOscillatorComplex.cc
 * 
 * \todo Provide optional key printing
 */
template<typename StateType, typename Derivs, typename ODE_Engine>
class Simulated
{
public:
  using StreamedArray=StateType;
  
  template <typename STATE>
  Simulated(STATE&& stateInit, Derivs derivs, std::initializer_list<std::string> keyLabels, ODE_Engine ode)
    : state_{std::forward<STATE>(stateInit)},
      derivs_{derivs},
      keyPrinter_{"Simulated",keyLabels},
      ode_{ode} {
        auto& labels{keyPrinter_.getLabels()};
        if (labels.size()<SubscriptLimit<StateType>::_(state_)) labels.insert(labels.end(),SubscriptLimit<StateType>::_(state_)-labels.size(),"N/A");
        else if (labels.size()>SubscriptLimit<StateType>::_(state_)) throw std::runtime_error("More keys than values in Simulated");
      }
  
  auto getTime() const {return t_;}

  void step(double deltaT, std::ostream& logStream) {ode_.step(deltaT,logStream,derivs_,t_,state_);}

  auto getDtDid() const {return ode_.getDtDid();}
  
  std::ostream& streamParameters(std::ostream& os) const {return ode_.streamParameters(os<<"\nSimulated.\n");}

  std::tuple<std::ostream&,StateType> stream(std::ostream& os, int precision) const
  {
    using namespace cppqedutils;
    for (size_t i=0; i < SubscriptLimit<StateType>::_(state_); i++) os<<FormDouble(precision)(Subscript_c<StateType>::_(state_,i))<<' ';
    return {os,state_.copy()};
  }

  iarchive& readFromArrayOnlyArchive(iarchive& iar) {return iar & state_;}

  /** structure of Simulated archives:
   * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry )
   */
  template <typename Archive>
  Archive& stateIO(Archive& ar) {return ode_.stateIO(ar & state_ & t_);} // state should precede time in order to be compatible with array-only archives
  
  std::ostream& streamKey(std::ostream& os) const {size_t i=3; return keyPrinter_.stream(os,i);}

  std::ostream& logOnEnd(std::ostream& os) const {return ode_.logOnEnd(os);}
  
private:
  double t_=0.;
  
  StateType state_;
  
  Derivs derivs_;
  
  KeyPrinter keyPrinter_;
  
  ODE_Engine ode_;
  
};


/// Deduction guide:
template<typename StateType, typename Derivs, typename ODE_Engine>
Simulated(StateType, Derivs, std::initializer_list<std::string>, ODE_Engine) -> Simulated<StateType,Derivs,ODE_Engine> ;



template <typename StateType, typename Derivs, typename ODE_Engine>
struct trajectory::MakeSerializationMetadata<Simulated<StateType,Derivs,ODE_Engine>>
{
  static auto _() {return SerializationMetadata{typeid(ElementType_t<StateType>).name(),"Simulated",Rank_v<StateType>};}
};


namespace simulated {

template<typename ODE_Engine, typename StateType, typename Derivs, typename ... ODE_EngineCtorParams>
auto make(StateType&& stateInit, Derivs derivs, std::initializer_list<std::string> keyLabels, ODE_EngineCtorParams&&... odePack)
{
  return Simulated<std::decay_t<StateType>,Derivs,ODE_Engine>{std::forward<StateType>(stateInit),derivs,keyLabels,ODE_Engine{std::forward<ODE_EngineCtorParams>(odePack)...}};
}

template<typename StateType, typename Derivs, typename ... ODE_EngineCtorParams>
auto makeBoost(StateType&& stateInit, Derivs derivs, std::initializer_list<std::string> keyLabels, ODE_EngineCtorParams&&... odePack)
{
  return make<ODE_EngineBoost<std::decay_t<StateType>>>(std::forward<StateType>(stateInit),derivs,keyLabels,std::forward<ODE_EngineCtorParams>(odePack)...);
}

} // simulated


} // cppqedutils


#endif // CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED
