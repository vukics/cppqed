// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED

#include "Master.h"

#include "DO_Display.tcc"
#include "QuantumTrajectory.h"

#include "MathExtensions.h"

#include "BlitzArraySliceIterator.tcc"
#include "ComplexArrayExtensions.h"

#include <boost/range/algorithm_ext/for_each.hpp>

namespace quantumtrajectory {


namespace master {


template<int RANK>
Base<RANK>::Base(DO_Ptr rho,
                 typename QuantumSystem::Ptr qs,
                 const master::Pars& p,
                 const DensityOperatorLow& scaleAbs
                 )
  : QuantumTrajectory(qs,true,
                      rho->getArray(),
                      bind(&Base<RANK>::derivs,this,_1,_2,_3),
                      initialTimeStep<RANK>(qs),
                      p.logLevel,p.epsRel,p.epsAbs,
                      scaleAbs,
                      evolved::MakerGSL<DensityOperatorLow>(p.sf,p.nextDtTryCorrectionFactor)),
    rho_(rho)
{
  if (!this->applicableInMaster()) throw master::SystemNotApplicable();
  QuantumTrajectory::checkDimension(rho);

}


template<int RANK>
void Base<RANK>::derivs(double t, const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow) const
{
  drhodtLow=0;

  binaryIter(rhoLow,drhodtLow,bind(&QuantumTrajectory::QuantumSystemWrapper::addContribution,this,t,_1,_2,this->getT0()));

  {
    linalg::CMatrix drhodtMatrixView(blitzplusplus::binaryArray(drhodtLow));
    linalg::calculateTwoTimesRealPartOfSelf(drhodtMatrixView);
  }

  // Now act with the reset operator --- implement this in terms of the individual jumps by iteration and addition

  for (size_t i=0; i<this->template nAvr<structure::LA_Li>(); i++) {
    try {
      this->getLi()->actWithSuperoperator(t,rhoLow,drhodtLow,i);
    } catch (const structure::SuperoperatorNotImplementedException&) {
      DensityOperatorLow rhotemp(rhoLow.copy());
      UnaryFunction functionLi(bind(&Liouvillean::actWithJ,this->getLi(),t,_1,i));
      unaryIter(rhotemp,functionLi);
      blitzplusplus::hermitianConjugateSelf(rhotemp);
      unaryIter(rhotemp,functionLi);
      drhodtLow+=rhotemp;
    }
  }

}

template<int RANK>
void 
Base<RANK>::step_v(double deltaT)
{
  this->getEvolved()->step(deltaT);
  if (const auto ex=this->getEx()) {
    using namespace blitzplusplus;
    DensityOperatorLow rhoLow(rho_->getArray());
    UnaryFunction functionEx(bind(&Exact::actWithU,ex,this->getTime(),_1,this->getT0()));
    unaryIter(rhoLow,functionEx);
    hermitianConjugateSelf(rhoLow);
    unaryIter(rhoLow,functionEx);
    rhoLow=conj(rhoLow);
  }

  QuantumTrajectory::setT0();


  // The following "smoothing" of rho_ has proven to be necessary for the algorithm to remain stable:
  // We make the approximately Hermitian and normalized rho_ exactly so.

  {
    linalg::CMatrix m(rho_->matrixView());
    linalg::calculateTwoTimesRealPartOfSelf(m); 
    // here we get two times of what is desired, but it is anyway renormalized in the next step
  }

  rho_->renorm();

}


template<int RANK>
std::ostream& Base<RANK>::displayParameters_v(std::ostream& os) const
{
  using namespace std;

  this->displayCharacteristics( this->getQS()->displayParameters( Adaptive::displayParameters_v(os)<<"Solving Master equation."<<addToParameterDisplay()<<endl<<endl ) )<<endl;

  if (const auto li=this->getLi()) {
    os<<"Decay channels:\n";
    {
      size_t i=0;
      li->displayKey(os,i);
    }
    os<<"Explicit superoperator calculations: ";
    DensityOperator rhotemp(rho_->getDimensions());
    {
      int n=0;
      for (int i=0; i<li->nAvr(); ++i)
        try {
          li->actWithSuperoperator(0,rho_->getArray(),rhotemp.getArray(),i);
          os<<i<<' '; ++n;
        }
        catch (const structure::SuperoperatorNotImplementedException&) {}
      if (!n) os<<"none";
    }
    os<<endl;

  }
  
  return os;
}


// NEEDS_WORK needs to create again the iterators because the array stored in the iterator IS NOT transposed because transpose does not touch the data. For this reason one must be very careful with reusing SliceIterators because probably most of the times it does not produce the required semantics. Alternatively, one could imagine a rebind functionality in SliceIterators.

template<int RANK>
void Base<RANK>::unaryIter(DensityOperatorLow& rhoLow, UnaryFunction function) const
{
  using namespace blitzplusplus::vfmsi;
  boost::for_each(fullRange<Left>(rhoLow),function);
}


template<int RANK>
void Base<RANK>::binaryIter(const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow, BinaryFunction function) const
{
  using namespace blitzplusplus::vfmsi;
  for_each(fullRange<Left>(rhoLow),fullRange<Left>(drhodtLow),function);
}


} // master


} // quantumtrajectory


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED
