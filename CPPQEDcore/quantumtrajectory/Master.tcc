// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED

#include "Master.h"

#include "QuantumTrajectory.h"

#include "MathExtensions.h"

#include "BlitzArraySliceIterator.tcc"
#include "ComplexArrayExtensions.h"

// #define USE_BOOST_PROGRESS_TIMER_PROFILING
#include "Profiling.h"

#include <boost/range/algorithm_ext/for_each.hpp>

namespace quantumtrajectory {


namespace master {


template<int RANK>
Base<RANK>::Base(DensityOperator& rho,
                 typename QuantumSystem::Ptr qs,
                 const master::Pars& p,
                 const DensityOperatorLow& scaleAbs
                 )
  : QuantumTrajectory(qs,true,
             rho.getArray(),
             bind(&Base<RANK>::derivs,this,_1,_2,_3),
             initialTimeStep<RANK>(qs),
             p,
             scaleAbs,
             evolved::MakerGSL<DensityOperatorLow>(p.sf,p.nextDtTryCorrectionFactor)),
    rho_(rho)
{
  if (!getQSW().applicableInMaster()) throw master::SystemNotApplicable();
  // If the interaction picture is non-unitary, the density matrix in IP is non-Hermitian. This cannot be allowed here because then the calculation of the Hamiltonian part of the dynamics as implemented below would fail.
  // if (!li_) throw master::NoLiouvillean();

  if (rho!=*getQSW().getQS()) throw DimensionalityMismatchException();

}


template<int RANK>
void Base<RANK>::derivs(double t, const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow) const
{
  drhodtLow=0;

  PROGRESS_TIMER_IN_POINT(getOstream());

  binaryIter(rhoLow,drhodtLow,bind(&QuantumTrajectory::QuantumSystemWrapper::addContribution,getQSW(),t,_1,_2,this->tIntPic0_));

  PROGRESS_TIMER_OUT_POINT("Hamiltonian");

  {
    PROGRESS_TIMER_IN_POINT(getOstream())
    linalg::CMatrix drhodtMatrixView(blitzplusplus::binaryArray(drhodtLow));
    linalg::calculateTwoTimesRealPartOfSelf(drhodtMatrixView);
    PROGRESS_TIMER_OUT_POINT("RealPartOfSelf")
  }

  // Now act with the reset operator --- implement this in terms of
  // the individual jumps by iteration and addition

  for (size_t i=0; i<getQSW().template nAvr<structure::LA_Li>(); i++) {
    PROGRESS_TIMER_IN_POINT( getOstream() )
    DensityOperatorLow rhotemp(rhoLow.copy());
    UnaryFunction functionLi(bind(&Liouvillean::actWithJ,getQSW().getLi(),t,_1,i));
    unaryIter(rhotemp,functionLi);
    blitzplusplus::hermitianConjugateSelf(rhotemp);
    unaryIter(rhotemp,functionLi);
    drhodtLow+=rhotemp;
    PROGRESS_TIMER_OUT_POINT("Liouvillean")
  }

}

template<int RANK>
void 
Base<RANK>::step_v(double deltaT)
{
  PROGRESS_TIMER_IN_POINT( getOstream() )
  getEvolved()->step(deltaT);
  PROGRESS_TIMER_OUT_POINT("Evolved step total")
    ;
  if (const auto ex=getQSW().getEx()) {
    PROGRESS_TIMER_IN_POINT( getOstream() )
    using namespace blitzplusplus;
    DensityOperatorLow rhoLow(rho_.getArray());
    UnaryFunction functionEx(bind(&Exact::actWithU,ex,getDtDid(),_1,this->tIntPic0_));
    unaryIter(rhoLow,functionEx);
    // rhoLow=hermitianConjugate(rhoLow) 
    hermitianConjugateSelf(rhoLow);
    unaryIter(rhoLow,functionEx);
    rhoLow=conj(rhoLow);
    PROGRESS_TIMER_OUT_POINT("Exact")
  }

  this->tIntPic0_=getTime();


  // The following "smoothing" of rho_ has proven to be necessary for the algorithm to remain stable:
  // We make the approximately Hermitian and normalized rho_ exactly so.

  PROGRESS_TIMER_IN_POINT( getOstream() ) ;

  {
    linalg::CMatrix m(rho_.matrixView());
    linalg::calculateTwoTimesRealPartOfSelf(m); 
    // here we get two times of what is desired, but it is anyway renormalized in the next step
  }

  rho_.renorm();

  PROGRESS_TIMER_OUT_POINT("Smoothing")

}


template<int RANK>
std::ostream& Base<RANK>::displayParameters_v(std::ostream& os) const
{
  using namespace std;

  getQSW().displayCharacteristics( getQSW().getQS()->displayParameters( Adaptive::displayParameters_v(os)<<"# Solving Master equation."<<addToParameterDisplay()<<endl<<endl ) )<<endl;

  if (const auto li=getQSW().getLi()) {
    os<<"# Decay channels:\n";
    {
      size_t i=0;
      li->displayKey(os,i);
    }
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


template<int RANK>
void BaseFast<RANK>::unaryIter(DensityOperatorLow& rhoLow, UnaryFunction function) const
{
  for_each(blitzplusplus::basi_fast::fullRange(rhoLow,slicesData_),function);
}


template<int RANK>
void BaseFast<RANK>::binaryIter(const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow, BinaryFunction function) const
{
  for_each(blitzplusplus::basi_fast::fullRange(   rhoLow,slicesData_),
           blitzplusplus::basi_fast::fullRange(drhodtLow,slicesData_),function);
}

} // master


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED
