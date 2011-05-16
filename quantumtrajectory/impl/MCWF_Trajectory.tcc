// -*- C++ -*-
#ifndef   _MCWF_TRAJECTORY_IMPL_H
#define   _MCWF_TRAJECTORY_IMPL_H

#include "ParsMCWF_Trajectory.h"

#include "StateVector.h"
#include "Structure.h"

#include "FormDouble.h"

#include <fstream>
#include <vector>

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
  }
}

/*
const vector<HS_Vector> ReadBasis(const string& basisFile, size_t basisDim, size_t systemDim)
{
  vector<HS_Vector> basis(basisDim);
  if (basisDim) {
    ifstream FILE(basisFile.data());
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
	 evolved::MakerGSL<StateVectorLow>(),
	 randomized::MakerGSL()),
    tIntPic0_(0),// random_(0),
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
    dpLimit_(p.dpLimit),
    svdc_(p.svdc),
    file_(p.ofn),
    initFile_(p.initFile+".sv")
{
  using namespace std;

  if (psi!=sys) throw DimensionalityMismatchException();
  

  if (initFile_!=".sv") {
    ifstream file(initFile_.data());
    if (!file.is_open()) throw MCWF_TrajectoryFileOpeningException(initFile_);
#define READ_INTO_PSI {StateVectorLow psiTemp; file>>psiTemp; psi_=psiTemp; psi_.renorm();}
    READ_INTO_PSI;
  }

  if (file_!="") {
    {
      ifstream file(file_.data());
      if (!file.is_open() || file.peek()==EOF) return;
    }
    {
      ifstream file((file_+".sv").data());
      if (!file.is_open()) throw MCWF_TrajectoryFileOpeningException(file_+".sv");
      {
	READ_INTO_PSI;
      }
#undef READ_INTO_PSI
      {
	char c;
	file>>c; // eat '#'
	if (c!='#') return;
      }
      double t0, dtTry;
      file>>t0; file>>dtTry;
      getEvolved()->setTime(t0); getEvolved()->setDtTry(dtTry);
    }

  }

}


template<int RANK>
MCWF_Trajectory<RANK>::~MCWF_Trajectory()
{
  using namespace std;

  if (file_!="") {
    ofstream file((file_+".sv").data());
    file<<psi_();
    file<<"# "<<getEvolved()->getTime()<<' '<<getEvolved()->getDtTry()<<endl;
  }

}


template<int RANK>
void MCWF_Trajectory<RANK>::displayMore(int precision) const
{
  using namespace std;

  ostream& os=getOstream();

  Averaged::display(getTime(),psi_,os,precision,av_);

  displayEvenMore(precision);

  os<</*'\t'<<high()(random_)<<*/endl;

  static unsigned svdCount_=0;
  if (svdc_ && !(svdCount_%svdc_)) os<<psi_();
  svdCount_++;
}


template<int RANK>
void MCWF_Trajectory<RANK>::step(double Dt) const
{
  // Coherent time development
  if (ha_) getEvolved()->step(Dt);
  else getEvolved()->update(getTime()+Dt,Dt);

  double t=getTime();

  if (ex_) {
    ex_->actWithU(getDtDid(),psi_());
    tIntPic0_=t;
  }

  psi_.renorm();

  // Jump
  if (li_) {

    typename Liouvillean::Probabilities probas(Liouvillean::probabilities(t,psi_,li_));

    std::vector<indexSVL_tuple> probasSpecial;

    {
      for (int i=0; i<probas.size(); i++)
	if (probas(i)<0) {
	  StateVector psiTemp(psi_);
	  Liouvillean::actWithJ(t,psiTemp(),i,li_);
	  probasSpecial.push_back(indexSVL_tuple(i,psiTemp()));

	  probas(i)=mathutils::sqr(psiTemp.renorm());

	  // psiTemp disappears here, but its storage does not, because the ownership is taken over by the SVL in the tuple

	}
    }
    

    double dpOverDt=std::accumulate(probas.begin(),probas.end(),0.);

    double random_=(*getRandomized())()/getDtDid();

    {
      int jumpNo=0;
      for (; random_>0 && jumpNo!=probas.size(); random_-=probas(jumpNo++))
	;

      if(random_<0) { // Jump No. jumpNo-1 occurs
	struct helper
	{
	  static bool p(int i, indexSVL_tuple j) {return i==j.template get<0>();} // NEEDS_WORK how to express this with lambda?
	};

	typename std::vector<indexSVL_tuple>::const_iterator i(find_if(probasSpecial,bind(&helper::p,--jumpNo,_1))); // See whether it's a special jump
	if (i!=probasSpecial.end())
	  // special jump
	  psi_()=i->template get<1>(); // RHS already normalized above
	else {
	  // normal  jump
	  Liouvillean::actWithJ(t,psi_(),jumpNo,li_); 
	  psi_()/=sqrt(probas(jumpNo));
	}
      }
    }

    if (dpOverDt*getEvolved()->getDtTry()>dpLimit_)
      getEvolved()->setDtTry(dpLimit_/dpOverDt);
    // Time-step management --- jump probability should not overshoot dpLimit

  }

}


template<int RANK>
void MCWF_Trajectory<RANK>::displayParameters() const
{
  using namespace std;
  Base::displayParameters();

  ostream& os=getOstream();

  os<<"# MCWF Trajectory Parameters: dpLimit="<<dpLimit_<<endl;
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
    typename Liouvillean::Probabilities probas(Liouvillean::probabilities(0,psi_,li_));
    
    for (int i=0; i<probas.size(); i++) if (probas(i)<0) os<<i<<' ';
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
