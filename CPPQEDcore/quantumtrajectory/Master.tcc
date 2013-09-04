// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED

#include "Master.h"

#include "Structure.h"

#include "MathExtensions.h"

#include "Algorithm.h"
#include "BlitzArraySliceIterator.tcc"
#include "ComplexArrayExtensions.h"
#include "Range.h"

// #define USE_BOOST_PROGRESS_TIMER_PROFILING
#include "Profiling.h"


namespace quantumtrajectory {


namespace master {


template<int RANK>
Base<RANK>::Base(DensityOperator& rho,
                 typename QuantumSystem::Ptr qs,
                 const master::Pars& p,
                 const DensityOperatorLow& scaleAbs
                 )
  : Adaptive(rho(),
             bind(&Base<RANK>::derivs,this,_1,_2,_3),
             trajectory::initialTimeStep(qs->highestFrequency()),
             p,
             scaleAbs,
             evolved::MakerGSL<DensityOperatorLow>(p.sf,p.nextDtTryCorretionFactor)),
    rho_(rho),
    tIntPic0_(0),
    qs_(qs)
{
  if (!qs_.isUnitary()) throw master::NonUnitaryIP();
  // If the interaction picture is non-unitary, the density matrix in IP is non-Hermitian. This cannot be allowed here because then the calculation of the Hamiltonian part of the dynamics as implemented below would fail.
  // if (!li_) throw master::NoLiouvillean();

  if (rho!=*qs_.getQS()) throw DimensionalityMismatchException();

}


template<int RANK>
void Base<RANK>::derivs(double t, const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow) const
{
  drhodtLow=0;

  PROGRESS_TIMER_IN_POINT(getOstream());

  binaryIter(rhoLow,drhodtLow,bind(&QuantumSystemWrapper::addContribution,qs_,t,_1,_2,tIntPic0_));

  PROGRESS_TIMER_OUT_POINT("Hamiltonian");

  {
    PROGRESS_TIMER_IN_POINT(getOstream())
    linalg::CMatrix drhodtMatrixView(blitzplusplus::binaryArray(drhodtLow));
    linalg::calculateTwoTimesRealPartOfSelf(drhodtMatrixView);
    PROGRESS_TIMER_OUT_POINT("RealPartOfSelf")
  }

  // Now act with the reset operator --- implement this in terms of
  // the individual jumps by iteration and addition

  for (size_t i=0; i<qs_.template nAvr<structure::LA_Li>(); i++) {
    PROGRESS_TIMER_IN_POINT( getOstream() )
    DensityOperatorLow rhotemp(rhoLow.copy());
    UnaryFunction functionLi(bind(&Liouvillean::actWithJ,qs_.getLi(),t,_1,i));
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
  if (const typename Exact::Ptr ex=qs_.getEx()) {
    PROGRESS_TIMER_IN_POINT( getOstream() )
    using namespace blitzplusplus;
    DensityOperatorLow rhoLow(rho_());
    UnaryFunction functionEx(bind(&Exact::actWithU,ex,getDtDid(),_1,tIntPic0_));
    unaryIter(rhoLow,functionEx);
    // rhoLow=hermitianConjugate(rhoLow) 
    hermitianConjugateSelf(rhoLow);
    unaryIter(rhoLow,functionEx);
    rhoLow=conj(rhoLow);
    PROGRESS_TIMER_OUT_POINT("Exact")
  }

  tIntPic0_=getTime();


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
  
  qs_.displayCharacteristics( qs_.getQS()->displayParameters( Adaptive::displayParameters_v(os)<<"# Solving Master equation."<<addToParameterDisplay()<<endl<<endl ) )<<endl;

  if (const typename Liouvillean::Ptr li=qs_.getLi()) {
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
  cpputils::for_each(fullRange<Left>(rhoLow),begin<Left>(drhodtLow),function);
}


template<int RANK>
void BaseFast<RANK>::unaryIter(DensityOperatorLow& rhoLow, UnaryFunction function) const
{
  boost::for_each(blitzplusplus::basi_fast::fullRange(rhoLow,slicesData_),function);
}


template<int RANK>
void BaseFast<RANK>::binaryIter(const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow, BinaryFunction function) const
{
  cpputils::for_each(blitzplusplus::basi_fast::fullRange(   rhoLow,slicesData_),
                     blitzplusplus::basi_fast::begin    (drhodtLow,slicesData_),function);
}

} // master

template<int RANK, typename V, bool IS_FAST>
cpputils::oarchive& Master<RANK,V,IS_FAST>::writeMeta_v(cpputils::oarchive& oar) const
{
  trajectory::SerializationMetadata meta;
  meta.rank=2*RANK;
  meta.trajectoryType="Master";
  return oar & meta;
}

} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_MASTER_TCC_INCLUDED
