#define BLITZ_ARRAY_LARGEST_RANK 15

#include "EvolutionComposite.h"

#include "impl/AveragingUtils.tcc"
#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"
#include "WignerFunction.h"


using namespace std;

typedef quantumdata::StateVector<5> StateVector;
typedef quantumtrajectory::TimeAveragingMCWF_Trajectory<5> Trajectory;

int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  qbit::ParsPumpedLossy pplqb(p);
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);
  
  double
    &wfLimitXL=p.add("wfLimitXL","",2.),
    &wfLimitXU=p.add("wfLimitXU","",2.),
    &wfLimitYL=p.add("wfLimitYL","",2.),
    &wfLimitYU=p.add("wfLimitYU","",2.),
    &wfStep=p.add("wfStep","",1.);
    
  int &wfCutoff=p.add("wfCutoff","",10);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  const qbit::Ptr qbit(qbit::make(pplqb,qmp));
  const mode::Ptr mode(mode::make<mode::AveragedQuadratures>(pplm ,qmp));
  const jaynescummings::Ptr jc(jaynescummings::make(qbit,mode,pjc)),
                            jcRDO(jaynescummings::make<ReducedDensityOperatorNegativity<2,tmptools::Vector<0> > >(qbit,mode,pjc,"JaynesCummings0-4",jc->getDimensions())),
                            jcCorr(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc));
  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*mode::init(pplm));

  Trajectory
    traj(psi,composite::make(Act<0,4>(jcCorr),Act<1,4>(jcRDO),Act<2,4>(jc),Act<3,4>(jc)),pe,pe.relaxationTime);
    
  {
    const string stateFileName(pe.ofn+".state");
    ifstream stateFile(stateFileName.c_str(),ios_base::binary);
    cpputils::iarchive stateArchive(stateFile);
    traj.readState(stateArchive);
  }
  
  quantumdata::DensityOperator<2> rho2(jcRDO->getDimensions());
  inflate(traj.getAverages()(blitz::Range(29,traj.getAverages().size()-2)),rho2,true);
 
  quantumdata::DensityOperator<1> rho(quantumdata::partialTrace<tmptools::Vector<1>,quantumdata::DensityOperator<1> >(rho2,quantumdata::densityOperatorize<1>));
  
  for (double x=wfLimitXL; x<wfLimitXU; x+=wfStep) {
    for (double y=wfLimitYL; y<wfLimitYU; y+=wfStep) cout<<x<<"\t"<<y<<"\t"<<cpputils::wignerFunction(rho,x,y,wfCutoff)<<endl;
    cout<<endl;
  }  

}


/*
typedef TTD_IDXTINY(2) Idx;
cerr<<averages(5000)<<" "<<averages(5001)<<" "<<rho(Idx(0,5),Idx(1,121))<<" "<<rho(Idx(1,121),Idx(0,5))<<endl;
cerr<<rho(37,29)<<" "<<rho2(Idx(0,37),Idx(0,29))+rho2(Idx(1,37),Idx(1,29))<<endl;
*/
