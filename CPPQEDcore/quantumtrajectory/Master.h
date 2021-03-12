// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED

#include "Structure.h"

#include "StreamDensityOperator.h"
#include "QuantumTrajectory.h"
#include "Types.h"

#include "ComplexArrayExtensions.h"
#include "MathExtensions.h"
#include "VectorFromMatrixSliceIterator.h"



namespace quantumtrajectory {


/// Auxiliary tools to Master
namespace master {


typedef trajectory::ParsEvolved Pars;


/// Thrown if the system is not applicable in Master-equation evolution
/**
 * \see structure::ExactCommon::applicableInMaster, \ref masterequationlimitations
 */
struct SystemNotApplicable : std::runtime_error {SystemNotApplicable() : std::runtime_error("") {}};


} // master



/// An \link trajectory::Adaptive Adaptive\endlink trajectory class representing Master equation evolution from a \link quantumdata::DensityOperator density-operator\endlink initial condition
/**
 * \see \ref masterequation.

 * \note The ODE driver underlying this class needs to store several (typically 6–7, this depends on the chosen driver) density-operator instants.
 *
 * \tparam RANK arity of the Hilbert space
 * \tparam V has the same function as the template parameter `V` in stream_densityoperator::_, which class is used here for deriving quantum averages to stream from the evolved density operator
 * 
 */
template<int RANK, typename ODE_Engine, typename V=tmptools::V_Empty>
class Master : public structure::QuantumSystemWrapper<RANK,true>
{
protected:
  Master(Master&&) = default; Master& operator=(Master&&) = default;

public:
  typedef structure::QuantumSystem<RANK> QuantumSystem;

  typedef typename quantumdata::Types<RANK>::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  Master(DensityOperator&& rho, ///< the density operator to be evolved
         typename QuantumSystem::Ptr sys, ///< object representing the quantum system
         ODE_Engine ode,
         bool negativity ///< governs whether entanglement should be calculated, cf. stream_densityoperator::_, quantumdata::negPT
         ) 
  : structure::QuantumSystemWrapper<RANK,true>{qs,true},
    rho_{std::forward<DensityOperator>(rho)},
    ode_(ode),
    dos_(this->getAv(),negativity)
  {
    if (!this->applicableInMaster()) throw master::SystemNotApplicable();
    QTraj::checkDimension(rho);
  }
    
  auto getTime() const {return t_;}

  void step(double deltaT, std::ostream& logStream);

  auto getDtDid() const {return ode_.getDtDid();}
  
  std::ostream& streamParameters(std::ostream& os) const;
  
  auto stream(std::ostream& os, int precision) const {return dos_.stream(t_,rho_,os,precision);}
  
  iarchive& readFromArrayOnlyArchive(iarchive& iar) {return rho_.readFromArrayOnlyArchive(iar);}

  /** structure of Master archives:
  * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry ) - t0
  */
  template <typename Archive>
  Archive& stateIO(Archive& ar) {return ode_.stateIO(ar & rho_ & t_) & t0_;} // state should precede time in order to be compatible with array-only archives

  std::ostream& streamKey_v(std::ostream& os) const {size_t i=3; return dos_.streamKey(os,i);}
  
  std::ostream& logOnEnd(std::ostream& os) const {return ode_.logOnEnd(os);}

private:
  double t_=0., t0_=0.;
  
  DensityOperator rho_;
  
  ODE_Engine ode_;

  const stream_densityoperator::_<RANK,V> dos_;

};


} // quantumtrajectory


template <int RANK, typename V>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::Master<RANK,V>>
{
  static auto _() {return SerializationMetadata{typeid(dcomp).name(),"Master",RANK};}
};



template<int RANK, typename ODE_Engine, typename V>
void quantumtrajectory::Master<RANK,V>::step(double deltaT, std::ostream& logStream)
{
  const auto derivs = [this](const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow, double t)
    {
      using namespace blitzplusplus::vfmsi;

      drhodtLow=0;
      
      for (auto&& [rhoS,drhodtS] : boost::combine(fullRange<Left>(rhoLow),fullRange<Left>(drhodtLow)) )
        this->addContribution(t,rhoS,drhodtS,this->getT0());
      
      {
        linalg::CMatrix drhodtMatrixView(blitzplusplus::binaryArray(drhodtLow));
        linalg::calculateTwoTimesRealPartOfSelf(drhodtMatrixView);
      }
      
      // Now act with the reset operator --- implement this in terms of the individual jumps by iteration and addition
      
      for (size_t i=0; i < this->template nAvr<structure::LA_Li>(); i++) {
        try {
          this->getLi()->actWithSuperoperator(t,rhoLow,drhodtLow,i);
        } catch (const structure::SuperoperatorNotImplementedException&) {
          DensityOperatorLow rhotemp(rhoLow.copy());
#define UNARY_ITER for (auto& rhoS : fullRange<Left>(rhotemp)) this->getLi()->actWithJ(t,rhoS,i);
          UNARY_ITER; blitzplusplus::hermitianConjugateSelf(rhotemp); UNARY_ITER;
#undef UNARY_ITER
          drhodtLow+=rhotemp;
        }
      }            
    };
  
  ode_.step(deltaT,logStream,derivs,t_,rho_.getArray());

  if (const auto ex=this->getEx()) {
    using namespace blitzplusplus;
    DensityOperatorLow rhoLow(rho_.getArray());
    for (auto& rhoS : vfmsi::fullRange<Left>(rhoLow)) ex->actWithU(this->getTime(),rhoS,this->getT0());
    hermitianConjugateSelf(rhoLow);
    for (auto& rhoS : vfmsi::fullRange<Left>(rhoLow)) ex->actWithU(this->getTime(),rhoS,this->getT0());
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
std::ostream& quantumtrajectory::Master<RANK,V>::streamParameters(std::ostream& os) const
{
  using namespace std;

  this->streamCharacteristics( this->getQS()->streamParameters(ode_.streamParameters(os)<<"Solving Master equation."<<endl<<endl ) )<<endl;

  if (const auto li=this->getLi()) {
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
        catch (const structure::SuperoperatorNotImplementedException&) {}
      if (!n) os<<"none";
    }
    os<<endl;

  }
  
  return os;
}



#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
