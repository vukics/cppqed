// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED

#include "Master.h"

#include "QuantumTrajectory.h"

#include "MathExtensions.h"

#include "VectorFromMatrixSliceIterator.h"
#include "ComplexArrayExtensions.h"


template<int RANK, typename V>
quantumtrajectory::Master<RANK,V>::Master(DensityOperator&& rho,
                                          typename QuantumSystem::Ptr qs,
                                          const master::Pars& p,
                                          bool negativity,
                                          const DensityOperatorLow& scaleAbs)
  : QTraj(qs,true,
          rho.getArray(),
          [this](double t, const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow)
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
          },
          initialTimeStep<RANK>(qs),
          p.logLevel,p.epsRel,p.epsAbs,
          scaleAbs,
          evolved::MakerGSL<DensityOperatorLow>(p.sf,p.nextDtTryCorrectionFactor)),
    rho_(std::move(rho)),
    dos_(this->getAv(),negativity)
{
  if (!this->applicableInMaster()) throw master::SystemNotApplicable();
  QTraj::checkDimension(rho);
}


template<int RANK, typename V>
void quantumtrajectory::Master<RANK,V>::step_v(double deltaT, std::ostream& logStream)
{
  this->getEvolved()->step(deltaT,logStream);
  if (const auto ex=this->getEx()) {
    using namespace blitzplusplus;
    using namespace vfmsi;
    DensityOperatorLow rhoLow(rho_.getArray());
    for (auto& rhoS : fullRange<Left>(rhoLow)) ex->actWithU(this->getTime(),rhoS,this->getT0());
    hermitianConjugateSelf(rhoLow);
    for (auto& rhoS : fullRange<Left>(rhoLow)) ex->actWithU(this->getTime(),rhoS,this->getT0());
    rhoLow=conj(rhoLow);
  }

  QTraj::setT0();


  // The following "smoothing" of rho_ has proven to be necessary for the algorithm to remain stable:
  // We make the approximately Hermitian and normalized rho_ exactly so.

  {
    linalg::CMatrix m(rho_.matrixView());
    linalg::calculateTwoTimesRealPartOfSelf(m); 
    // here we get two times of what is needed, but it is anyway renormalized in the next step
  }

  rho_.renorm();

}


template<int RANK, typename V>
std::ostream& quantumtrajectory::Master<RANK,V>::streamParameters_v(std::ostream& os) const
{
  using namespace std;

  this->streamCharacteristics( this->getQS()->streamParameters(Adaptive::streamParameters_v(os)<<"Solving Master equation."<<endl<<endl ) )<<endl;

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


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED
