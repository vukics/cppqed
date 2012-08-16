#include "Evolution.h"

#include "benchmarking.h"

#include "Mode_.h"

using namespace mode;


int main(int, char**)
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsPumpedLossy pplm(p);

  pplm.eta=dcomp(4.,-2.);
  pplm.cutoff=1000000;

  // ****** ****** ****** ****** ****** ******

  const PumpedLossyModeSch<false> mode(pplm);

  // benchmark(mode,mode,Vector<0>());

  mode.displayParameters(std::cout);

  quantumdata::StateVector<1> psi(mode.getDimensions());

  {
    Randomized::SmartPtr Ran(MakerGSL()(1001));
    boost::generate(psi(),bind(&Randomized::dcompRan,Ran));
  }

  // 1

  // if (doDisplay) cout<<blitzplusplus::unaryArray(psi());
  quantumdata::StateVector<1> psiout(psi);

  {
    const structure::Hamiltonian<1>* ha=dynamic_cast<const structure::Hamiltonian<1>*>(&mode);

    boost::progress_timer t;
    for (unsigned count=0; count<1000; ++count)
      structure::Hamiltonian<1>::addContribution(0.,psi(),psiout(),0.,ha);
  }

  // if (doDisplay) cout<<blitzplusplus::unaryArray(psiout());


}


