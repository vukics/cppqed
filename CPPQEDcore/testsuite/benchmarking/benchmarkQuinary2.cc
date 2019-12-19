// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "benchmarking.h"

#include "Mode_.h"

using namespace mode;


struct Ia : structure::Interaction<2>, structure::TridiagonalHamiltonian<2,false>
{
  Ia(const ModeBase* m) : structure::Interaction<2>(Frees(m,m)), structure::TridiagonalHamiltonian<2,false>(mode::aop(m)*mode::aop(m).dagger()) {}
};


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  const size_t cutoff=16;

  // ****** ****** ****** ****** ****** ******

  const ModeBase mode(cutoff);

  const Ia ia(&mode);
  const structure::Interaction<5> dummy(structure::Interaction<5>::Frees(&mode,&mode,&mode,&mode,&mode));

  const Composite<composite::result_of::make_vector<Act<0,2>,Act<0,1,2,3,4> >::type>
    sys(makeComposite(Act<0,2>(ia),
                      Act<0,1,2,3,4>(dummy)
                                     ));

  benchmark(sys,ia,tmptools::Vector<0,2>());

}


