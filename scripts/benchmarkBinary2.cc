#include "EvolutionHigh.h"

#include "benchmarking.h"

#include "Mode.h"

using namespace mode;


struct Ia : structure::Interaction<2>, structure::TridiagonalHamiltonian<2,false>
{
  Ia(const ModeBase* m) : structure::Interaction<2>(Frees(m,m)), structure::TridiagonalHamiltonian<2,false>(mode::aop(m)*mode::aop(m).dagger()) {}
};


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  const size_t cutoff=1000;

  // ****** ****** ****** ****** ****** ******


  const ModeBase mode(cutoff);

  const Ia ia(&mode);

  const Composite<composite::result_of::make_vector<Act<0,1> >::type> sys(makeComposite(Act<0,1>(ia)));

  benchmark(sys);

}


