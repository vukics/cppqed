#include "Composite.h"

#include "Randomized.h"

#include <boost/progress.hpp>


using namespace std;
using namespace randomized;


template<int RANK>
void benchmark(const structure::QuantumSystem<RANK>& sys, bool doDisplay=false, size_t nRepeat=1000)
{
  const structure::Hamiltonian<RANK>* ha=dynamic_cast<const structure::Hamiltonian<RANK>*>(&sys);

  sys.displayParameters(std::cout);

  quantumdata::StateVector<RANK> psi(sys.getDimensions());

  {
    Randomized::SmartPtr Ran(MakerGSL()(1001));
    boost::generate(psi(),bind(&Randomized::dcompRan,Ran));
  }

  if (doDisplay) cout<<blitzplusplus::rankOneArray(psi());

  quantumdata::StateVector<RANK> psiout(psi);

  {
    boost::progress_timer t;
    for (unsigned count=0; count<nRepeat; ++count)
      structure::Hamiltonian<RANK>::addContribution(0.,psi(),psiout(),0.,ha);
  }

  if (doDisplay) cout<<blitzplusplus::rankOneArray(psiout());

}
