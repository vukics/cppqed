// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_IMPL_MCWF_TRAJECTORY_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_IMPL_MCWF_TRAJECTORY_TCC_INCLUDED

#include "MCWF_Trajectory.h"

#include "ParsMCWF_Trajectory.h"

#include "StateVector.h"
#include "impl/StochasticTrajectory.tcc"
#include "Structure.h"

#include "impl/FormDouble.tcc"
#include "SmartPtr.h"

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/complex.hpp>
#endif // DO_NOT_USE_BOOST_SERIALIZATION
#include <fstream>


// #include<blitz/tinyvec-et.h>

// #include "Eigen.h"




namespace quantumtrajectory {

using namespace trajectory;

//////////////////
//
// Implementations
//
//////////////////

template<int RANK>
void MCWF_Trajectory<RANK>::derivs(double t, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  if (const typename Hamiltonian::Ptr ha=qs_.getHa()) {
    dpsidt=0;

    ha->addContribution(t,psi,dpsidt,tIntPic0_);
    logger_.hamiltonianCalled();
  }
}

/*
const vector<HS_Vector> ReadBasis(const string& basisFile, size_t basisDim, size_t systemDim)
{
  vector<HS_Vector> basis(basisDim);
  if (basisDim) {
    ifstream FILE(basisFile.c_str());
    for (vector<HS_Vector>::iterator i=basis.begin(); i!=basis.end(); ++i) {
      HS_Vector temp(systemDim);
      FILE>>temp;
      (*i)=temp;
    }
  }
  return basis;
}

*/

template<int RANK> template<typename SYS>
MCWF_Trajectory<RANK>::MCWF_Trajectory(
				       StateVector& psi,
				       const SYS& sys,
				       const ParsMCWF_Trajectory& p,
				       const StateVectorLow& scaleAbs
				       )
  : TrajectoryBase(p),
    Base(psi(),
	 bind(&MCWF_Trajectory::derivs,this,_1,_2,_3),
	 1./(cpputils::sharedPointerize(sys)->highestFrequency()*Base::factor()),
	 scaleAbs,
	 p,
	 evolved::MakerGSL<StateVectorLow>(p.sf,p.nextDtTryCorretionFactor),
	 randomized::MakerGSL()),
    tIntPic0_(0),
    psi_(psi),
    qs_(cpputils::sharedPointerize(sys),Base::noise()),
    dpLimit_(p.dpLimit), overshootTolerance_(p.overshootTolerance),
    svdc_(p.svdc),
    firstSVDisplay_(p.firstSVDisplay),
    svdPrecision_(p.svdPrecision ? p.svdPrecision : getPrecision()),
    svdCount_(0),
#ifndef DO_NOT_USE_BOOST_SERIALIZATION
    binarySVFile_(p.binarySVFile),
    svExtension_(binarySVFile_?".svbin":".sv"),
#else // DO_NOT_USE_BOOST_SERIALIZATION
    svExtension_(".sv"),
#endif // DO_NOT_USE_BOOST_SERIALIZATION
    file_(p.ofn),
    initFile_(p.initFile+svExtension_),
    logger_(p.logLevel,qs_.getHa(),getOstream())
{
  using namespace std;

  if (psi!=*qs_.getQS()) throw DimensionalityMismatchException();
  
  if (initFile_!=svExtension_) {
    ifstream file(initFile_.c_str());
    if (!file.is_open()) throw MCWF_TrajectoryFileOpeningException(initFile_);
    readState(file,true);
  }

  class NoPresetDtTry {};

  try {

  if (file_!="") {
    {
      ifstream file(file_.c_str());
      if (!file.is_open() || file.peek()==EOF) throw NoPresetDtTry();
    }
    {
      ifstream file((file_+svExtension_).c_str());
      if (!file.is_open()) throw MCWF_TrajectoryFileOpeningException(file_+svExtension_);
      readState(file);
    }
  }
  else
    throw NoPresetDtTry();

  } catch (NoPresetDtTry) {
    if (p.logLevel>1) getOstream()<<"# Adjusting initial dtTry\n";
    manageTimeStep(qs_.probabilities(0.,psi_),getEvolved().get(),false);
    // Initially, dpLimit should not be overshot, either.
  }
}


template<int RANK>
MCWF_Trajectory<RANK>::~MCWF_Trajectory()
{
  using namespace std;

  if (file_!="") {
    ofstream file((file_+svExtension_).c_str());
    writeState(file);
  }

}


template<int RANK>
void MCWF_Trajectory<RANK>::readState(std::ifstream &ifs, bool onlySV)
{
  using namespace std;
  
  StateVectorLow psiTemp;
  double t0, dtTry;
#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  if (!binarySVFile_) {
#endif // DO_NOT_USE_BOOST_SERIALIZATION
    ifs>>psiTemp;
    psi_=psiTemp;
    psi_.renorm();
    if (onlySV) return;
    
#define EAT_COMMENT_CHAR  {char c; ifs>>c;  if (c!='#') throw MCWF_TrajectoryFileParsingException(file_+svExtension_);} ifs.exceptions ( ifstream::failbit | ifstream::badbit | ifstream::eofbit ); \
/**/
    EAT_COMMENT_CHAR
    ifs>>*getRandomized();
    EAT_COMMENT_CHAR
#undef EAT_COMMENT_CHAR
    ifs>>t0>>dtTry;
#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  }
  else {
    boost::archive::binary_iarchive ia(ifs);
    for(int i=0, d; i<=RANK; i++) ia>>d; // eat RANK and dimensions
    ia>>psiTemp;
    psi_=psiTemp;
    if (onlySV) return;
    ia>>t0>>dtTry;
    ifs>>*getRandomized();
  }
#endif // DO_NOT_USE_BOOST_SERIALIZATION
  getEvolved()->update(t0,dtTry); getEvolved()->setDtDid(0); svdCount_=1;
  if (qs_.getEx()) tIntPic0_=t0;
  getOstream()<<"# Next timestep to try: "<<dtTry<<endl;
}


template<int RANK>
void MCWF_Trajectory<RANK>::writeState(std::ofstream &ofs) const
{
  using namespace std;

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  if (!binarySVFile_) {
#endif // DO_NOT_USE_BOOST_SERIALIZATION
    ofs<<formdouble::zeroWidth(getPrecision())(psi_())<<"\n# "
       <<*getRandomized()<<"\n# "
       <<formdouble::positive(getPrecision())(getTime ())
       <<formdouble::positive(getPrecision())(getDtTry())<<endl;
#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  }
  else {
    boost::archive::binary_oarchive oa(ofs);
    double t=getTime(); double dttry=getDtTry();
    int r=RANK;
    oa<<r;
    for(int i=0;i<RANK;i++) {
      int d=psi_().extent(i);
      oa<<d;
    }
    oa<<psi_()<<t<<dttry;
    ofs<<*getRandomized();
  }
#endif // DO_NOT_USE_BOOST_SERIALIZATION

}



template<int RANK>
void MCWF_Trajectory<RANK>::displayMore() const
{
  using namespace std;

  ostream& os=getOstream();

  qs_.display(getTime(),psi_,os,getPrecision());

  displayEvenMore();

  os<<endl;

  if (svdc_ && !(svdCount_%svdc_) && (svdCount_||firstSVDisplay_)) os<<FormDouble(svdPrecision_,0)(psi_());
  svdCount_++;
}


// NEEDS_WORK factor out the functions coherentTimeDevelopment calculateDpOverDtSpecialSet manageTimeStep performJump

template<int RANK>
double MCWF_Trajectory<RANK>::coherentTimeDevelopment(double Dt) const
{
  if (qs_.getHa()) getEvolved()->step(Dt);
  else {
    double stepToDo=qs_.getLi() ? std::min(getDtTry(),Dt) : Dt;
    getEvolved()->update(getTime()+stepToDo,getDtTry());
    logger_.logFailedSteps(getEvolved()->nFailedSteps());
  }

  double t=getTime();

  if (const typename Exact::Ptr ex=qs_.getEx()) {
    ex->actWithU(getDtDid(),psi_());
    tIntPic0_=t;
  }

  logger_.processNorm(psi_.renorm());

  return t;
}


template<int RANK>
const typename MCWF_Trajectory<RANK>::IndexSVL_tuples
MCWF_Trajectory<RANK>::calculateDpOverDtSpecialSet(DpOverDtSet* dpOverDtSet, double t) const
{
  IndexSVL_tuples res;
  for (int i=0; i<dpOverDtSet->size(); i++)
    if ((*dpOverDtSet)(i)<0) {
      StateVector psiTemp(psi_);
      qs_.actWithJ(t,psiTemp(),i);
      res.push_back(IndexSVL_tuple(i,psiTemp()));
      (*dpOverDtSet)(i)=mathutils::sqr(psiTemp.renorm());
    } // psiTemp disappears here, but its storage does not, because the ownership is taken over by the SVL in the tuple
  return res;
}


template<int RANK>
bool MCWF_Trajectory<RANK>::manageTimeStep(const DpOverDtSet& dpOverDtSet, evolved::TimeStepBookkeeper* evolvedCache, bool logControl) const
{
  const double dpOverDt=std::accumulate(dpOverDtSet.begin(),dpOverDtSet.end(),0.);
  const double dtDid=getDtDid(), dtTry=getDtTry();

  // Assumption: overshootTolerance_>=1 (equality is the limiting case of no tolerance)
  if (dpOverDt*dtDid>overshootTolerance_*dpLimit_) {
    evolvedCache->setDtTry(dpLimit_/dpOverDt);
    (*getEvolved())=*evolvedCache;
    logger_.stepBack(dpOverDt*dtDid,dtDid,getDtTry(),getTime(),logControl);
    return true; // Step-back required.
  }
  else if (dpOverDt*dtTry>dpLimit_) {
    // dtTry-adjustment for next step required
    getEvolved()->setDtTry(dpLimit_/dpOverDt);
    logger_.overshot(dpOverDt*dtTry,dtTry,getDtTry(),logControl);
  }
  
  return false; // Step-back not required.
}


template<int RANK>
void MCWF_Trajectory<RANK>::performJump(const DpOverDtSet& dpOverDtSet, const IndexSVL_tuples& dpOverDtSpecialSet, double t) const
{
  double random=(*getRandomized())()/getDtDid();

  int jumpNo=0;
  for (; random>0 && jumpNo!=dpOverDtSet.size(); random-=dpOverDtSet(jumpNo++))
    ;

  if(random<0) { // Jump No. jumpNo-1 occurs
    struct helper
    {
      static bool p(int i, IndexSVL_tuple j) {return i==j.template get<0>();} // NEEDS_WORK how to express this with lambda?
    };

    typename IndexSVL_tuples::const_iterator i(find_if(dpOverDtSpecialSet,bind(&helper::p,--jumpNo,_1))); // See whether it's a special jump
    if (i!=dpOverDtSpecialSet.end())
      // special jump
      psi_()=i->template get<1>(); // RHS already normalized above
    else {
      // normal  jump
      qs_.actWithJ(t,psi_(),jumpNo);
      double normFactor=sqrt(dpOverDtSet(jumpNo));
      if (!boost::math::isfinite(normFactor)) throw structure::InfiniteDetectedException();
      psi_()/=normFactor;
    }
	
    logger_.jumpOccured(t,jumpNo);
  }
}


template<int RANK>
void MCWF_Trajectory<RANK>::step(double Dt) const
{
  const StateVectorLow psiCache(psi_().copy());
  evolved::TimeStepBookkeeper evolvedCache(*getEvolved()); // This cannot be const since dtTry might change.

  double t=coherentTimeDevelopment(Dt);

  if (const typename Liouvillean::Ptr li=qs_.getLi()) {

    DpOverDtSet dpOverDtSet(li->probabilities(t,psi_));
    IndexSVL_tuples dpOverDtSpecialSet=calculateDpOverDtSpecialSet(&dpOverDtSet,t);

    while (manageTimeStep(dpOverDtSet,&evolvedCache)) {
      psi_()=psiCache;
      t=coherentTimeDevelopment(Dt); // the next try
      dpOverDtSet=li->probabilities(t,psi_);
      dpOverDtSpecialSet=calculateDpOverDtSpecialSet(&dpOverDtSet,t);
    }

    // Jump
    performJump(dpOverDtSet,dpOverDtSpecialSet,t);

  }

  logger_.step();

}


template<int RANK>
void MCWF_Trajectory<RANK>::displayParameters() const
{
  using namespace std;
  Base::displayParameters();

  ostream& os=getOstream();

  os<<"# MCWF Trajectory Parameters: dpLimit="<<dpLimit_<<" (overshoot tolerance factor)="<<overshootTolerance_<<endl;
  if (svdc_) os<<"# Displaying State Vector in every "<<svdc_<<" Display"<<endl;
  os<<endl;

  qs_.getQS()->displayParameters(os);

  os<<"# System characteristics: "
    <<(qs_.getEx() ? "Interaction picture, "   : "")
    <<(qs_.getHa() ? "Hamiltonian evolution, " : "")
    <<(qs_.getLi() ? "Liouvillean evolution, " : "")
    <<(qs_.getAv() ? "calculates Averages."    : "")
    <<endl;

  if (const typename Liouvillean::Ptr li=qs_.getLi()) {
    os<<"# Decay channels:\n";
    {
      size_t i=0;
      li->displayKey(getOstream(),i);
    }
    os<<"# Alternative jumps: ";
    {
      const DpOverDtSet dpOverDtSet(li->probabilities(0,psi_));
      int n=0;
      for (int i=0; i<dpOverDtSet.size(); i++) if (dpOverDtSet(i)<0) {os<<i<<' '; n++;}
      if (!n) os<<"none";
    }
    os<<endl;
  }

}


template<int RANK>
size_t MCWF_Trajectory<RANK>::displayMoreKey() const
{
  size_t i=3;
  qs_.displayAveragedKey(getOstream(),i);
  return i;
}


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_IMPL_MCWF_TRAJECTORY_TCC_INCLUDED
