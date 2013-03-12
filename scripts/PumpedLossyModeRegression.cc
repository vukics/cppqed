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

  string& referenceStateFileName=p.add("referenceStateFileName","",string());
  
  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  pe.nTraj=10000//00
    ;
  const list<size_t> sections=boost::assign::list_of(100)(200)(500)(1000)(2000)(5000)(10000)//(20000)(50000)(100000)
    ;

  Ptr modeMaster(make(pplm,QMP_UIP)), mode(make(pplm,QMP_IP));
  StateVector psi(init(pplm));

  DensityOperator rho(psi);

  Master      <1> m(rho,*modeMaster,pe,false);

  EnsembleMCWF<1> e(psi,*mode      ,pe,false);
  e.displayParameters(cout);

  ifstream referenceStateFile(referenceStateFileName.c_str());
  
  while ( (referenceStateFile.peek(), !referenceStateFile.eof()) ) {
    {
      { // restoring next state of m
        string buffer; streamsize n; referenceStateFile>>n; buffer.resize(n);
        referenceStateFile.read(&buffer[0],n);
        istringstream iss(buffer,ios_base::binary);
        cpputils::iarchive referenceStateArchive(iss);
        m.readState(referenceStateArchive);
      }
      cout<<m.getTime()<<'\t';
      // evolving e to the same instant:
      e.evolve(m.getTime()-e.getTime());
    }
    for (list<size_t>::const_iterator i=sections.begin(); i!=sections.end(); ++i) {
      double avr=0;
      for (size_t begin=0; begin<pe.nTraj; begin+=*i)
        avr+=frobeniusNorm(rho-e.averageInRange(begin,*i));
      cout<<avr*(double(*i)/pe.nTraj)<<'\t';
    }
    cout<<endl;
  }




}



// PumpedLossyModeConvergence --nTh 10 --cutoff 100 --kappa 5 --eta "(-10,9)" --T 0.5 --Dt 0.0005 --precision 6 --minitFock 20 --seed 1000 --logLevel 1 --logLevel -1
