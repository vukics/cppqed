#include "EvolutionHigh.h"

#include "benchmarking.h"

#include "Mode.h"

using namespace mode;
using namespace blitzplusplus;

struct Ia : structure::Interaction<10>, structure::Hamiltonian<10,structure::NO_TIME>
{
  Ia(const ModeBase* m) : structure::Interaction<10>(Frees(m,m,m,m,m,m,m,m,m,m)),
			  h_(aop(m)*aop(m).dagger()*aop(m)*aop(m).dagger()*aop(m)*aop(m).dagger()) {}
  
  

  void addContribution(const StateVectorLow& psi, StateVectorLow& dpsidt) const
  {
    cpputils::for_each(basi::fullRange(psi,v),basi::begin(dpsidt,v),bind(&quantumoperator::apply<6>,_1,_2,h_));
  }

private:
  static const tmptools::Vector<0,2,4,5,6,8> v;

  const quantumoperator::Tridiagonal<6> h_;
};


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  const size_t cutoff=4;

  // ****** ****** ****** ****** ****** ******

  const ModeBase mode(cutoff);

  const Ia ia(&mode);

  const Composite<composite::result_of::make_vector<Act<0,1,2,3,4,5,6,7,8,9> >::type>
    sys(makeComposite(Act<0,1,2,3,4,5,6,7,8,9>(ia)));

  benchmark(sys);

}


