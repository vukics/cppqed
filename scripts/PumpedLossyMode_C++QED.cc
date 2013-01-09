#include "Evolution.h"
#include "Mode.h"

#include "BichromaticMode.h"
#include "AveragingUtils.h"


using namespace std ;
using namespace mode;


typedef averagingUtils::Collecting<1> Collecting;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  ParsBichromatic pplm(p); 

  bool
    &alternative=p.add("alternative","Alternative mode",false),
    &doDisplay=p.add("doDisplay","Display diagonal elements of density operator",false);

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if (pe.evol==EM_MASTER && qmp==QMP_IP) qmp=QMP_UIP;

  Collecting::Collection collection; collection.push_back(new AveragedQuadratures());
  if (doDisplay) collection.push_back(new ReducedDensityOperator<1>("",pplm.cutoff));

  Ptr mode(alternative ? Ptr(new PumpedLossyModeIP_NoExact(pplm)) : make<Collecting>(pplm,qmp,collection));

  StateVector psi(mode::init(pplm));

  evolve(psi,mode,pe);




}

