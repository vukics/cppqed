#include "EvolutionHigh.h"

#include "benchmarking.h"

#include "Mode.h"

using namespace mode;


struct Ia : structure::Interaction<2>, structure::TridiagonalHamiltonian<2,false>
{
  Ia(const ModeBase* m) : structure::Interaction<2>(Frees(m,m)), structure::TridiagonalHamiltonian<2,false>(mode::aop(m)*mode::aop(m).dagger()) {}
};


int main(int, char**)
{
  // ****** Parameters of the Problem

  const size_t cutoff=4;

  // ****** ****** ****** ****** ****** ******

  const ModeBase mode(cutoff);

  const Ia ia(&mode);
  const structure::Interaction<10> dummy(structure::Interaction<10>::Frees(&mode,&mode,&mode,&mode,&mode,&mode,&mode,&mode,&mode,&mode));

  const Composite<composite::result_of::make_vector<Act<0,2>,Act<0,1,2,3,4,5,6,7,8,9> >::type>
    sys(makeComposite(Act<0,2>(ia),
		      Act<0,1,2,3,4,5,6,7,8,9>(dummy)
				     ));

  benchmark(sys,ia,Vector<0,2>());

}


