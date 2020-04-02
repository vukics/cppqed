// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "benchmarking.h"

#include "Mode_.h"

using namespace mode;


struct Ia : structure::Interaction<2>, structure::TridiagonalHamiltonian<2,false>
{
  Ia(const ModeBase* m) : structure::Interaction<2>(Frees(m,m)), structure::TridiagonalHamiltonian<2,false>(mode::aop(m)*mode::aop(m).dagger()) {}
};


int main(int, char*[])
{
  // ****** Parameters of the Problem

  const size_t cutoff=1000;

  // ****** ****** ****** ****** ****** ******


  const ModeBase mode(cutoff);

  const Ia ia(&mode);

  const Composite<composite::result_of::make_vector<Act<0,1> >::type> sys(makeComposite(Act<0,1>(ia)));

  benchmark(sys,ia,Vector<0,1>());

}


