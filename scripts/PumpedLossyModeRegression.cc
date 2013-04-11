#include "Evolution.h"
#include "Mode.h"


using namespace std ;
using namespace mode;



int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  ParsPumpedLossy pplm(p); 

  string
    &referenceStateFileName=p.add("referenceStateFileName","",string()), 
    &         stateFileName=p.add("stateFileName"         ,"",string());

  size_t& nSections=p.add("nSections","",size_t(10));
  
  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  const vector<size_t> sections=boost::assign::list_of(100)(200)(500)(1000)(2000)(5000)(10000)(20000)(50000)(100000)
    ;
  pe.nTraj=sections[nSections-1];
    ;

  Ptr modeMaster(make(pplm,QMP_UIP)), mode(make(pplm,QMP_IP));
  StateVector psi(init(pplm));

  DensityOperator rho(psi);

  Master      <1> m(rho,*modeMaster,pe,false);

  EnsembleMCWF<1> e(psi,*mode      ,pe,false);

  ifstream referenceStateFile(referenceStateFileName.c_str());

  { // Determine whether there is already a saved previous state
    ifstream stateFile(stateFileName.c_str());
    
    if ( stateFile.is_open() && (stateFile.peek(), !stateFile.eof()) ) {
      trajectory::readViaSStream(e,stateFile,true);
      // Restoring m’s state to the first time instant greater than e’s time:
      while ( (referenceStateFile.peek(), !referenceStateFile.eof()) && m.getTime()<e.getTime() ) trajectory::readViaSStream(m,referenceStateFile,false);
      cout<<"# Continuing up to time "<<pe.T<<endl;
    }
    else {
      e.displayParameters(cout);
      cout<<"\n# Running up to time "<<pe.T<<endl<<endl;
    }
  }

  while ( (referenceStateFile.peek(), !referenceStateFile.eof()) && e.getTime()<pe.T ) {
    {
      trajectory::readViaSStream(m,referenceStateFile,false);
      // evolving e to the same instant:
      e.evolve(m.getTime()-e.getTime());
      cout<<m.getTime()<<'\t'<<e.getDtDid()<<'\t';
    }
    for (vector<size_t>::const_iterator i=sections.begin(); i!=sections.begin()+nSections; ++i) {
      double avr=0;
      size_t nTrajCurrent=*i;
      for (size_t begin=0; begin<pe.nTraj; begin+=nTrajCurrent)
        avr+=frobeniusNorm(rho-e.averageInRange(begin,nTrajCurrent));
      cout<<avr*(double(nTrajCurrent)/pe.nTraj)<<'\t';
    }
    cout<<endl;
  }

  e.logOnEnd(cout);

  if (stateFileName!="") { // Save final state
    ofstream stateFile(stateFileName.c_str(),ios_base::app);
    trajectory::writeViaSStream(e,&stateFile);  
  }  

}



// PumpedLossyModeConvergence --nTh 10 --cutoff 100 --kappa 5 --eta "(-10,9)" --T 0.5 --Dt 0.0005 --precision 6 --minitFock 20 --seed 1000 --logLevel 1 --logLevel -1
