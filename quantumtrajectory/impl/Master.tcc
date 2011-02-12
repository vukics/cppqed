// -*- C++ -*-
#ifndef   _MASTER_IMPL_H
#define   _MASTER_IMPL_H

#include "Structure.h"

#include "MathExtensions.h"

#include "Algorithm.h"
#include "ComplexArrayExtensions.h"
#include "Range.h"

// #define USE_BOOST_PROGRESS_TIMER_PROFILING
#include "Profiling.h"


namespace quantumtrajectory {


namespace master {


template<int RANK>
Base<RANK>::Base(DensityOperator& rho,
		 const QuantumSystem& sys,
		 const master::Pars& pt,
		 const DensityOperatorLow& scaleAbs
		 )
  : trajectory::TrajectoryBase(pt),
    TrajectoryBase(rho(),
		   bind(&Base<RANK>::derivs,this,_1,_2,_3),
		   1./(sys.highestFrequency()*TrajectoryBase::factor()),
		   scaleAbs,
		   pt,
		   evolved::MakerGSL<DensityOperatorLow>()),
    rho_(rho),
    tIntPic0_(0),
    qs_(&sys),
    ex_(structure::qse(&sys)),
    ha_(structure::qsh(&sys)),
    li_(structure::qsl(&sys))
{
  if (!Exact::isUnitary(ex_)) throw master::NonUnitaryIP();
  // If the interaction picture is non-unitary, the density matrix in
  // IP is non-Hermitian. This cannot be allowed here because then the
  // calculation of the Hamiltonian part of the dynamics as
  // implemented below would fail.
  // if (!li_) throw master::NoLiouvillean();

  if (rho!=sys) throw DimensionalityMismatchException();

}


template<int RANK>
void Base<RANK>::derivs(double t, const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow) const
{
  drhodtLow=0;

  PROGRESS_TIMER_IN_POINT(getOstream());

  binaryIter(rhoLow,drhodtLow,bind(&Hamiltonian::addContribution,t,_1,_2,tIntPic0_,ha_,structure::theStaticOne));

  PROGRESS_TIMER_OUT_POINT("Hamiltonian");

  {
    PROGRESS_TIMER_IN_POINT(getOstream())
    linalg::CMatrix drhodtMatrixView(blitzplusplus::rankTwoArray(drhodtLow));
    linalg::calculateTwoTimesRealPartOfSelf(drhodtMatrixView);
    PROGRESS_TIMER_OUT_POINT("RealPartOfSelf")
  }

  // Now act with the reset operator --- implement this in terms of
  // the individual jumps by iteration and addition

  for (size_t i=0; i<Liouvillean::nJumps(li_); i++) {
    PROGRESS_TIMER_IN_POINT( getOstream() )
    DensityOperatorLow rhotemp(rhoLow.copy());
    UnaryFunction functionLi(bind(&Liouvillean::actWithJ,_1,i,li_,structure::theStaticOne));
    unaryIter(rhotemp,functionLi);
    blitzplusplus::hermitianConjugateSelf(rhotemp);
    unaryIter(rhotemp,functionLi);
    drhodtLow+=rhotemp;
    PROGRESS_TIMER_OUT_POINT("Liouvillean")
  }

}


template<int RANK>
void 
Base<RANK>::step(double deltaT) const
{
  PROGRESS_TIMER_IN_POINT( getOstream() )
  getEvolved()->step(deltaT);
  PROGRESS_TIMER_OUT_POINT("Evolved step total")

  if (ex_) {
    PROGRESS_TIMER_IN_POINT( getOstream() )
    using namespace blitzplusplus;
    DensityOperatorLow rhoLow(rho_());
    UnaryFunction functionEx(bind(&Exact::actWithU,getDtDid(),_1,ex_,structure::theStaticOne));
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
void
Base<RANK>::displayParameters() const
{
  using namespace std;
  TrajectoryBase::displayParameters();

  getOstream()<<"# Solving Master equation."<<addToParameterDisplay()<<endl<<endl;

  qs_->displayParameters(getOstream());

}


// NEEDS_WORK needs to create again the iterators because the
// array stored in the iterator IS NOT transposed because
// transpose does not touch the data. For this reason one must be
// very careful with reusing SliceIterators because probably most
// of the times it does not produce the required
// semantics. Alternatively, one could imagine a rebind
// functionality in SliceIterators.
template<int RANK>
void Base<RANK>::unaryIter(DensityOperatorLow& rhoLow, UnaryFunction function) const
{
  boost::for_each(fullRange(rhoLow,blitzplusplus::vfmsi::Left()),function);
}


template<int RANK>
void Base<RANK>::binaryIter(const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow, BinaryFunction function) const
{
  cpputils::for_each(fullRange(rhoLow,blitzplusplus::vfmsi::Left()),begin(drhodtLow,blitzplusplus::vfmsi::Left()),function);
}


template<int RANK>
void BaseFast<RANK>::unaryIter(DensityOperatorLow& rhoLow, UnaryFunction function) const
{
  using namespace blitzplusplus;
  boost::for_each(vfmsi_fast::fullRange(slicesData_,rhoLow,vfmsi::Left()),function);
}


template<int RANK>
void BaseFast<RANK>::binaryIter(const DensityOperatorLow& rhoLow, DensityOperatorLow& drhodtLow, BinaryFunction function) const
{
  using namespace blitzplusplus;
  cpputils::for_each(vfmsi_fast::fullRange(slicesData_,rhoLow,vfmsi::Left()),vfmsi_fast::begin(slicesData_,drhodtLow,vfmsi::Left()),function);
}



} // master

} // quantumtrajectory


#endif // _MASTER_IMPL_H
