// -*- C++ -*-
#ifndef   _MCWF_TRAJECTORY_IMPL_H
#define   _MCWF_TRAJECTORY_IMPL_H

#include "ParsMCWF_Trajectory.h"

#include "StateVector.h"
#include "Structure.h"

#include "FormDouble.h"

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
  if (ha_) {
    dpsidt=0;

    Hamiltonian::addContribution(t,psi,dpsidt,tIntPic0_,ha_);
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

template<int RANK>
MCWF_Trajectory<RANK>::MCWF_Trajectory(
				       StateVector& psi,
				       const QuantumSystem& sys,
				       const ParsMCWF_Trajectory& p,
				       const StateVectorLow& scaleAbs
				       )
  : TrajectoryBase(p),
    Base(psi(),
	 bind(&MCWF_Trajectory::derivs,this,_1,_2,_3),
	 1./(sys.highestFrequency()*Base::factor()),
	 scaleAbs,
	 p,
	 evolved::MakerGSL<StateVectorLow>(p.sf,p.nextDtTryCorretionFactor),
	 randomized::MakerGSL()),
    tIntPic0_(0),
    psi_(psi),
    qs_(&sys),
    ex_(structure::qse(&sys)),
    ha_(structure::qsh(&sys)),
    li_(Base::noise() 
	? 
	structure::qsl(&sys) 
	: 
	0),
    av_(structure::qsa(&sys)),
    dpLimit_(p.dpLimit), overshootTolerance_(p.overshootTolerance),
    svdc_(p.svdc),
    firstSVDisplay_(p.firstSVDisplay),
    svdCount_(0),
    file_(p.ofn),
    initFile_(p.initFile+".sv"),
    logger_(p.logLevel,ha_,getOstream())
{
  using namespace std;

  if (psi!=sys) throw DimensionalityMismatchException();
  
  if (initFile_!=".sv") {
    ifstream file(initFile_.c_str());
    if (!file.is_open()) throw MCWF_TrajectoryFileOpeningException(initFile_);
#define READ_INTO_PSI {StateVectorLow psiTemp; file>>psiTemp; psi_=psiTemp; psi_.renorm();}
    READ_INTO_PSI;
  }

  class NoPresetDtTry {};

  try {

  if (file_!="") {
    {
      ifstream file(file_.c_str());
      if (!file.is_open() || file.peek()==EOF) throw NoPresetDtTry();
    }
    {
      ifstream file((file_+".sv").c_str());
      if (!file.is_open()) throw MCWF_TrajectoryFileOpeningException(file_+".sv");

      {
	READ_INTO_PSI;
      }
#undef READ_INTO_PSI
      {
	char c;
	file>>c; // eat '#'
	if (c!='#') throw MCWF_TrajectoryFileParsingException(file_+".sv");
      }
      file.exceptions ( ifstream::failbit | ifstream::badbit | ifstream::eofbit );
      double t0, dtTry;
      file>>t0; file>>dtTry;
      getOstream()<<"# Next timestep to try: "<<dtTry<<std::endl;
      getEvolved()->update(t0,dtTry); getEvolved()->setDtDid(0); svdCount_=1;
      if (ex_) tIntPic0_=t0;
    }

  }
  else
    throw NoPresetDtTry();

  } catch (NoPresetDtTry) {
    if (p.logLevel>1) getOstream()<<"# Adjusting initial dtTry\n";
    manageTimeStep(Liouvillean::probabilities(0.,psi_,li_),getEvolved().get(),false);
    // Initially, dpLimit should not be overshot, either.
  }
}


template<int RANK>
MCWF_Trajectory<RANK>::~MCWF_Trajectory()
{
  using namespace std;

  if (file_!="") {
    ofstream file((file_+".sv").c_str());
    file<<psi_();
    file<<"\n# "<<getTime()<<' '<<getDtTry()<<endl;
  }

}


template<int RANK>
void MCWF_Trajectory<RANK>::displayMore(int precision) const
{
  using namespace std;

  ostream& os=getOstream();

  Averaged::display(getTime(),psi_,os,precision,av_);

  displayEvenMore(precision);

  os<<endl;

  if (svdc_ && !(svdCount_%svdc_) && (svdCount_||firstSVDisplay_)) os<<psi_();
  svdCount_++;
}



template<int RANK>
double MCWF_Trajectory<RANK>::coherentTimeDevelopment(double Dt) const
{
  if (ha_) getEvolved()->step(Dt);
  else {
    double stepToDo=li_ ? std::min(getDtTry(),Dt) : Dt;
    getEvolved()->update(getTime()+stepToDo,getDtTry());
    logger_.logFailedSteps(getEvolved()->nFailedSteps());
  }

  double t=getTime();

  if (ex_) {
    ex_->actWithU(getDtDid(),psi_());
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
      Liouvillean::actWithJ(t,psiTemp(),i,li_);
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
      Liouvillean::actWithJ(t,psi_(),jumpNo,li_);
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

  if (li_) {

    DpOverDtSet dpOverDtSet(Liouvillean::probabilities(t,psi_,li_));
    IndexSVL_tuples dpOverDtSpecialSet=calculateDpOverDtSpecialSet(&dpOverDtSet,t);

    while (manageTimeStep(dpOverDtSet,&evolvedCache)) {
      psi_()=psiCache;
      t=coherentTimeDevelopment(Dt); // the next try
      dpOverDtSet=Liouvillean::probabilities(t,psi_,li_);
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

  qs_->displayParameters(os);

  os<<"# System characteristics: "
    <<(ex_ ? "Interaction picture, "   : "")
    <<(ha_ ? "Hamiltonian evolution, " : "")
    <<(li_ ? "Liouvillean evolution, " : "")
    <<(av_ ? "calculates Averages."    : "")
    <<endl;


  os<<"# Alternative jumps: ";
  {
    const DpOverDtSet dpOverDtSet(Liouvillean::probabilities(0,psi_,li_));
    for (int i=0; i<dpOverDtSet.size(); i++) if (dpOverDtSet(i)<0) os<<i<<' ';
  }
  os<<endl;

}


template<int RANK>
size_t MCWF_Trajectory<RANK>::displayMoreKey() const
{
  size_t i=3;
  Averaged::displayKey(getOstream(),i,av_);
  return i;
}


} // quantumtrajectory


#endif // _MCWF_TRAJECTORY_IMPL_H
