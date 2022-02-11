// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "Particle.h"

using namespace std;
using namespace particle;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  evolution::Pars<> pe   (p); // Driver Parameters
  ParsPumped    ppart(p); 

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  Ptr part(make(ppart,qmp));

  if (!ppart.init.getSig() && !ppart.vClass) {cerr<<"Incorrect initial condition"<<endl; abort();}

  evolve(init(ppart),part,pe);


}

