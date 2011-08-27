#include "EvolutionHigh.h"

#include "Composite.h"

#include "Mode.h"
#include "Interaction.h"

#include "Randomized.h"

#include <boost/progress.hpp>

using namespace std ;
using namespace mode;
using namespace randomized;

const unsigned nRepeat=3000;


struct Ia : structure::Interaction<2>, structure::TridiagonalHamiltonian<2,false>
{
  Ia(const ModeBase* m) : structure::Interaction<2>(Frees(m,m)), structure::TridiagonalHamiltonian<2,false>(mode::aop(m)*mode::aop(m).dagger()) {}
};


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  const size_t cutoff=1000;

  //  update(p,argc,argv,"--");

  // ****** ****** ****** ****** ****** ******

  Randomized::SmartPtr Ran(MakerGSL()(1001));

  const ModeBase mode(cutoff);

  Ia ia(&mode);

  Composite<composite::result_of::make_vector<Act<0,1> >::type> sys(makeComposite(Act<0,1>(ia)));

  // const BinarySystem sys(ia);

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


