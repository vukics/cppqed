// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "benchmarking.h"

#include "Mode_.h"

using namespace mode;


struct Ia : structure::Interaction<3>, structure::TridiagonalHamiltonian<3,false>
{
  Ia(const ModeBase* m) : structure::Interaction<3>(Frees(m,m,m)), structure::TridiagonalHamiltonian<3,false>(mode::aop(m)*mode::aop(m).dagger()*mode::aop(m)) {}
};


int main(int, char*[])
{
  // ****** Parameters of the Problem

  const size_t cutoff=16;

  // ****** ****** ****** ****** ****** ******

  const ModeBase mode(cutoff);

  const Ia ia(&mode);
  const structure::Interaction<5> dummy(structure::Interaction<5>::Frees(&mode,&mode,&mode,&mode,&mode));

  const Composite<composite::result_of::make_vector<Act<0,2,3>,Act<0,1,2,3,4> >::type>
    sys(makeComposite(Act<0,2,3>(ia),
                      Act<0,1,2,3,4>(dummy)
                                     ));

  benchmark(sys,ia,tmptools::Vector<0,2,3>());

}


