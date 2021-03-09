// Copyright András Vukics 2021. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_UTILS_ODE_GSL_H_INCLUDED
#define CPPQEDCORE_UTILS_ODE_GSL_H_INCLUDED

#include "ArrayTraits.h"
#include "ODE.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include <cstddef>

namespace cppqedutils {
  
  
namespace ode_engine {

  
template <typename StateType>
class ControlledErrorStepperGSL
{
public:
  using state_type=StateType;
  using deriv_type=state_type;
  using time_type=double;
  using value_type=ElementType_t<StateType>;
  using stepper_category=bno::controlled_stepper_tag;
  
  ControlledErrorStepperGSL(double epsAbs, double epsRel)
    : step_{},
      con_{gsl_odeiv2_control_standard_new(epsAbs,epsRel,1,1)},
      yerr_{},
      dydt_out_{} {}
  
private:
  using Aux_t=std::tuple<SystemFunctional_t<double,StateType>&,Extents_t<Rank_v<StateType>>>;
  
public:
  /// Realizes part of the logic of [gsl_odeiv2_evolve_apply](http://www.gnu.org/software/gsl/doc/html/ode-initval.html#c.gsl_odeiv2_evolve_apply)
  /**
   * Differences:
   * - First-same-as-last steppers not supported (at the moment only RKCK is supported)
   * - There is no guard against the stepsize becoming infinitesimal 
   */
  template <typename System, typename S>
  auto try_step(System sys, S&& stateInOut, double& time, double& dtTry)
  {
    auto totalExtent=Size<StateType>::_(stateInOut);
    
    if (!yerr_.size()) { // this means that tryStep is called for the first time
      step_=gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck,totalExtent);
      yerr_.resize(totalExtent);
      dydt_out_.resize(totalExtent);
    }
#ifndef NDEBUG
    else if (!IsStorageContiguous<StateType>::_(stateInOut)) {
      throw NonContiguousStorageException{"ControlledErrorStepperGSL::tryStep"};
    }
    else if (yerr_.size()!=Size<StateType>::_(stateInOut) || dydt_out_.size()!=Size<StateType>::_(stateInOut)) {
      throw std::runtime_error("Dimensionality mismatch in ControlledErrorStepperGSL::tryStep");
    }
#endif // NDEBUG
    
    SystemFunctional_t<double,StateType> sysVariable=sys; // it’s important to convert `sys` (which can be a lambda) to a variable of known type
    
    Aux_t aux{sysVariable,Extents<StateType>::_(stateInOut)};
    
    gsl_odeiv2_system dydt{ControlledErrorStepperGSL::lowLevelSystemFunction,NULL,totalExtent,&aux};
    
    auto step_status = gsl_odeiv2_step_apply (step_, time, dtTry, Data<StateType>::_(stateInOut), yerr_.data(), NULL, dydt_out_.data(), &dydt);
    
    if (step_status == GSL_EFAULT || step_status == GSL_EBADFUNC) throw std::runtime_error("GSL bad step_status");

    if (step_status != GSL_SUCCESS) {
      dtTry *= 0.5;
      return bno::fail;
    }
    else {
      time+=dtTry;
      /* const auto hadjust_status = */
      gsl_odeiv2_control_hadjust (con_, step_, Data<StateType>::_(stateInOut), yerr_.data(), dydt_out_.data(), &dtTry);
      // if (hadjust_status == GSL_ODEIV_HADJ_DEC) throw std::runtime_error{"Unexpected change of dtTry in ControlledErrorStepperGSL::tryStep"};
      return bno::success;
    }
  }
  
private:
  static int lowLevelSystemFunction(double t, const double* y, double* dydt, void* aux)
  {
    auto [highLevelSystemFunction,arrayExtents] = *static_cast<Aux_t*>(aux);

    const StateType yInterface(Create_c<StateType>::_(y,arrayExtents));

    StateType dydtInterface(Create<StateType>::_(dydt,arrayExtents));

    highLevelSystemFunction(yInterface,dydtInterface,t);

    return GSL_SUCCESS;
  }
  
  gsl_odeiv2_step * step_;

  gsl_odeiv2_control * const con_;
  
  // make this movable easily, as ControlledErrorStepperGSL is carried around by value
  std::vector<value_type> yerr_, dydt_out_;
  
};


template <typename StateType>
struct MakeControlledErrorStepper<ControlledErrorStepperGSL<StateType>>
{
  static auto _(double epsRel, double epsAbs)
  {
    return std::make_tuple(ControlledErrorStepperGSL<StateType>{epsRel,epsAbs},
                           ControlledErrorStepperParameters<ControlledErrorStepperGSL<StateType>>{epsRel,epsAbs});
  }

  template <typename ParsBase>
  static auto _(const Pars<ParsBase>& p)
  {
    return _(p.epsRel,p.epsAbs);
  }

  
};


template <typename StateType>
inline const std::string StepperDescriptor<ControlledErrorStepperGSL<StateType>> = "GSL RKCK controlled stepper";


} // ode_engine


template <typename StateType>
using ODE_EngineGSL = ode_engine::Base<ode_engine::ControlledErrorStepperGSL<StateType>>;


} // cppqedutils

#endif // CPPQEDCORE_UTILS_ODE_GSL_H_INCLUDED
