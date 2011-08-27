#include "EvolutionHigh.h"

#include "Composite.h"

#include "Mode.h"
#include "Interaction.h"

#include "Randomized.h"

#include <boost/progress.hpp>

using namespace std ;
using namespace mode;
using namespace randomized;

const unsigned nRepeat=1000;


struct Ia : structure::Interaction<5>, structure::TridiagonalHamiltonian<5,false>
{
  Ia(const ModeBase* m) : structure::Interaction<5>(Frees(m,m,m,m,m)), structure::TridiagonalHamiltonian<5,false>(aop(m)*aop(m).dagger()*aop(m)*aop(m).dagger()*aop(m)) {}
};


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  const size_t cutoff=16;

  //  update(p,argc,argv,"--");

  // ****** ****** ****** ****** ****** ******

  Randomized::SmartPtr Ran(MakerGSL()(1001));

  const ModeBase mode(cutoff);

  Ia ia(&mode);

  Composite<composite::result_of::make_vector<Act<0,1,2,3,4> >::type>
    sys(makeComposite(Act<0,1,2,3,4>(ia)));

  static_cast<const structure::QuantumSystem<5>&>(sys).displayParameters(cout);

  quantumdata::StateVector<5> psi(sys.getDimensions());
  boost::generate(psi(),bind(&Randomized::dcompRan,Ran));

  // cout<<blitzplusplus::rankOneArray(psi());

  quantumdata::StateVector<5> psiout(psi);

  {
    boost::progress_timer t;
    for (unsigned count=0; count<nRepeat; ++count)
      structure::Hamiltonian<5>::addContribution(0.,psi(),psiout(),0.,&sys);
  }

  // cout<<blitzplusplus::rankOneArray(psiout());

}


