// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumTrajectory.h"


namespace quantumtrajectory {


/// Auxiliary tools to Master
namespace master {


typedef cppqedutils::ode::Pars<> Pars;


} // master



/// Master equation evolution from a \link quantumdata::DensityOperator density-operator\endlink initial condition
/**
 * \see \ref masterequation.
 * 
 * \note The ODE driver underlying this class needs to store several (typically 6–7, this depends on the chosen driver) density-operator instants.
 *
 */
template <size_t RANK,
          cppqedutils::ode::engine ODE_Engine,
          auto axesOfSubsystem = cppqedutils::retainedAxes<> // the axes belonging to one of the subsystems defined for entanglementMeasuresCalculation
          >
class Master
{
public:
  Master(Master&&) = default; Master& operator=(Master&&) = default;

  using StreamedArray=::structure::EV_Array;

  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  template <typename StateVector_OR_DensityOperator>
  Master(::structure::QuantumSystemPtr<RANK> sys, ///< object representing the quantum system
         StateVector_OR_DensityOperator&& state, ///< the state vector or density operator to be evolved
         ODE_Engine ode,
         EntanglementMeasuresSwitch ems ///< governs which entanglement measures should be calculated
         ) 
  : sys_{sys}, rho_{std::forward<StateVector_OR_DensityOperator>(state)}, ode_(ode), dos_(::structure::castAv(sys),ems)
  {
    if (!::structure::applicableInMaster(sys)) throw master::SystemNotApplicable();
    if (rho_!=*sys) throw DimensionalityMismatchException("during QuantumTrajectory construction");
  }

  auto getTime() const {return t_;}

  void step(double deltaT, std::ostream& logStream);

  auto getDtDid() const {return ode_.getDtDid();}
  
  std::ostream& streamParameters(std::ostream& os) const;
  
  auto stream(std::ostream& os, int precision) const {return dos_(t_,rho_,os,precision);}
  
  auto& readFromArrayOnlyArchive(cppqedutils::iarchive& iar) {DensityOperatorLow temp; iar & temp; rho_.getArray().reference(temp); return iar;}

  /** structure of Master archives:
  * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry )
  */
  // state should precede time in order to be compatible with array-only archives
  auto& stateIO(cppqedutils::iarchive& iar)
  {
    DensityOperatorLow temp;
    ode_.stateIO(iar & temp & t_);
    rho_.getArray().reference(temp);
    t0_=t_; // A very important step!
    return iar;
  }
  
  auto& stateIO(cppqedutils::oarchive& oar) {return ode_.stateIO(oar & rho_.getArray() & t_);}

  std::ostream& streamKey(std::ostream& os) const {size_t i=3; return dos_.streamKey(os,i);}
  
  std::ostream& logOnEnd(std::ostream& os) const {return ode_.logOnEnd(os);}

private:
  double t_=0., t0_=0.;
  
  ::structure::QuantumSystemPtr<RANK> sys_;
  
  DensityOperator rho_;
  
  ODE_Engine ode_;

  const DensityOperatorStreamer<RANK,V> dos_;

};


/// Deduction guides (note: `V` cannot be deduced this way, and partial deduction is not possible as of C++17):
template<typename System, int RANK, typename ODE_Engine>
Master(System, quantumdata::DensityOperator<RANK>, ODE_Engine, bool) -> Master<RANK,ODE_Engine,tmptools::V_Empty>;

template<typename System, int RANK, typename ODE_Engine>
Master(System, quantumdata::StateVector<RANK>, ODE_Engine, bool) -> Master<RANK,ODE_Engine,tmptools::V_Empty>;


namespace master {

template<typename ODE_Engine, typename V, typename StateVector_OR_DensityOperator>
auto make(::structure::QuantumSystemPtr<std::decay_t<StateVector_OR_DensityOperator>::N_RANK> sys,
          StateVector_OR_DensityOperator&& state, const Pars& p, EntanglementMeasuresSwitch ems)
{
  return Master<std::decay_t<StateVector_OR_DensityOperator>::N_RANK,ODE_Engine,V>(
    sys,std::forward<StateVector_OR_DensityOperator>(state),ODE_Engine{initialTimeStep(sys),p},ems);
}
  
} // master


} // quantumtrajectory


template <int RANK, typename ODE_Engine, typename V>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::Master<RANK,ODE_Engine,V>>
{
  static auto _() {return SerializationMetadata{"CArray","Master",2*RANK};}
};



template<int RANK, typename ODE_Engine, typename V>
void quantumtrajectory::Master<RANK,ODE_Engine,V>::step(double deltaT, std::ostream& logStream)
{
  // auto derivs = ;
  
  ode_.step(deltaT,logStream,[this](const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow, double t)
    {
      using namespace blitzplusplus::vfmsi;

      drhodtLow=0;
      
      for (auto&& [rhoS,drhodtS] : boost::combine(fullRange<Left>(rhoLow),fullRange<Left>(drhodtLow)) )
        ::structure::addContribution(sys_,t,rhoS,drhodtS,t0_);
      
      {
        linalg::CMatrix drhodtMatrixView(blitzplusplus::binaryArray(drhodtLow));
        linalg::calculateTwoTimesRealPartOfSelf(drhodtMatrixView);
      }
      
      // Now act with the reset operator --- implement this in terms of the individual jumps by iteration and addition
      
      for (size_t i=0; i < ::structure::nAvr<::structure::LA_Li>(sys_); i++) {
        try {
          ::structure::actWithSuperoperator(sys_,t,rhoLow,drhodtLow,i);
        } catch (const ::structure::SuperoperatorNotImplementedException&) {
          DensityOperatorLow rhotemp(rhoLow.copy());
#define UNARY_ITER for (auto& rhoS : fullRange<Left>(rhotemp)) ::structure::actWithJ(sys_,t,rhoS,i);
          UNARY_ITER; blitzplusplus::hermitianConjugateSelf(rhotemp); UNARY_ITER;
#undef UNARY_ITER
          drhodtLow+=rhotemp;
        }
      }            
    },t_,rho_.getArray());

  if (const auto ex=::structure::castEx(sys_)) {
    using namespace blitzplusplus; using namespace vfmsi;
    DensityOperatorLow rhoLow(rho_.getArray());
    for (auto& rhoS : fullRange<Left>(rhoLow)) ex->actWithU(this->getTime(),rhoS,t0_);
    hermitianConjugateSelf(rhoLow);
    for (auto& rhoS : fullRange<Left>(rhoLow)) ex->actWithU(this->getTime(),rhoS,t0_);
    rhoLow=conj(rhoLow);
  }

  t0_=t_;

  // The following "smoothing" of rho_ has proven necessary for the algorithm to remain stable:
  // We make the approximately Hermitian and normalized rho_ exactly so.

  {
    linalg::CMatrix m(rho_.matrixView());
    linalg::calculateTwoTimesRealPartOfSelf(m); 
    // here we get two times of what is needed, but it is anyway renormalized in the next step
  }

  rho_.renorm();

}


template<int RANK, typename ODE_Engine, typename V>
std::ostream& quantumtrajectory::Master<RANK,ODE_Engine,V>::streamParameters(std::ostream& os) const
{
  using namespace std;

  ::structure::streamCharacteristics(sys_, sys_->streamParameters(ode_.streamParameters(os)<<"Solving Master equation."<<endl<<endl ) )<<endl;

  if (const auto li=::structure::castLi(sys_)) {
    os<<"Decay channels:\n";
    {
      size_t i=0;
      li->streamKey(os,i);
    }
    os<<"Explicit superoperator calculations: ";
    DensityOperator rhotemp(rho_.getDimensions());
    {
      int n=0;
      for (int i=0; i<li->nAvr(); ++i)
        try {
          li->actWithSuperoperator(0,rho_.getArray(),rhotemp.getArray(),i);
          os<<i<<' '; ++n;
        }
        catch (const ::structure::SuperoperatorNotImplementedException&) {}
      if (!n) os<<"none";
    }
    os<<endl;

  }
  
  return os;
}



#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
