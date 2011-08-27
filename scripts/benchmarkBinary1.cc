#include "EvolutionHigh.h"

#include "BinarySystem.h"

#include "Mode.h"
#include "Interaction.h"

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
  pplm.cutoff=1000;

  //  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  Randomized::SmartPtr Ran(MakerGSL()(1001));

  PumpedLossyModeSch<false> modeH(pplm);
  ModeBase mode0(pplm.cutoff);

  const BinarySystem sys(structure::Interaction<2>(structure::Interaction<2>::Frees(&modeH,&mode0)));
  
  static_cast<const structure::QuantumSystem<2>&>(sys).displayParameters(cout);

  quantumdata::StateVector<2> psi(sys.getDimensions());
  boost::generate(psi(),bind(&Randomized::dcompRan,Ran));

  // cout<<psi();

  quantumdata::StateVector<2> psiout(psi);

  {
    boost::progress_timer t;
    for (unsigned count=0; count<nRepeat; ++count)
      structure::Hamiltonian<2>::addContribution(0.,psi(),psiout(),0.,&sys);
  }

  // cout<<psiout();

}


