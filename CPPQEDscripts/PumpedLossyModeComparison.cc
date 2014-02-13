#include "Evolution.h"
#include "Mode.h"


using namespace std ;
using namespace mode;


template<int RANK>
double compareDifferentSizes(const quantumdata::DensityOperator<RANK>& rho1, const quantumdata::DensityOperator<RANK>& rho2)
{
  if (all(rho1.getDimensions()==rho2.getDimensions())) return frobeniusNorm(rho1-rho2);
  else {
    const linalg::CMatrix m1(rho1.matrixView()), m2(rho2.matrixView());
    const blitz::Range common(0,min(m1.ubound(0),m2.ubound(0)));
    return frobeniusNorm(DensityOperator(DensityOperator::DensityOperatorLow(m1(common,common)-m2(common,common)),quantumdata::byReference));
  }
}


int main(int argc, char* argv[])
{
  struct TimeMismatch {TimeMismatch(double t1, double t2, double avr) {cerr<<t1<<" "<<t2<<" "<<avr<<endl;}};
  
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars pe(p); // Driver Parameters
  ParsPumpedLossy pplm1(p,"1"), pplm2(p,"2");

  string
    &stateFileName1=p.add("stateFile1","",string()), 
    &stateFileName2=p.add("stateFile2","",string());

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  Ptr modeMaster1(make(pplm1,QMP_UIP)), modeMaster2(make(pplm2,QMP_UIP));

  DensityOperator rho1(modeMaster1->getDimensions()), rho2(modeMaster2->getDimensions());

  Master<1> m1(rho1,*modeMaster1,pe,false), m2(rho2,*modeMaster2,pe,false);

  ifstream stateFile1(stateFileName1.c_str()), stateFile2(stateFileName2.c_str());

  double avr=0.; size_t c=0;
  while ( (stateFile1.peek(), !stateFile1.eof()) && (stateFile2.peek(), !stateFile2.eof()) ) {
    trajectory::readViaSStream(m1,stateFile1); trajectory::readViaSStream(m2,stateFile2);
    if (m1.getTime()!=m2.getTime()) throw TimeMismatch(m1.getTime(),m2.getTime(),avr/c);
    avr+=compareDifferentSizes(m1.getRho(),m2.getRho());
    c++;
  }
  cout<<avr/c<<endl;

}



// PumpedLossyModeConvergence --nTh 10 --cutoff 100 --kappa 5 --eta "(-10,9)" --T 0.5 --Dt 0.0005 --precision 6 --minitFock 20 --seed 1000 --logLevel 1 --logLevel -1
