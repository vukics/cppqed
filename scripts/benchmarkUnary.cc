#include "EvolutionHigh.h"

#include "Mode.h"

#include "Randomized.h"

#include <boost/progress.hpp>

using namespace std ;
using namespace mode;
using namespace randomized;

const unsigned nRepeat=1000;

int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsPumpedLossy pplm(p);

  pplm.eta=dcomp(4.,-2.);
  pplm.cutoff=1000000;

  //  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  Randomized::SmartPtr Ran(MakerGSL()(1001));

  PumpedLossyModeSch<false> mode(pplm);
  mode.displayParameters(std::cout);

  StateVector psi(mode.getTotalDimension());
  boost::generate(psi(),bind(&Randomized::dcompRan,Ran));

  // cout<<psi();

  StateVector psiout(psi);

  {
    boost::progress_timer t;
    for (unsigned count=0; count<nRepeat; ++count)
      structure::Hamiltonian<1>::addContribution(0.,psi(),psiout(),0.,&mode);
  }

  // cout<<psiout();

}


