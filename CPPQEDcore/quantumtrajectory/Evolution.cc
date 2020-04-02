// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution_.h"

#include "Pars.tcc"

#include<iostream>


std::ostream& evolution::operator<<(std::ostream& os, Method em)
{
  switch (em) {
    // case CONVERGENCE: return os<<"convergence";
  case SINGLE      : return os<<"single"      ;
  case ENSEMBLE    : return os<<"ensemble"    ;
  case MASTER      : return os<<"master"      ;
  case MASTER_FAST :        os<<"master_fast" ;
  }
  return os;
}


std::istream& evolution::operator>>(std::istream& is, Method& em) 
{
  Method emtemp=MASTER_FAST;
  std::string s;

  is>>s;
  /*if      (s=="convergence") emtemp=CONVERGENCE;
    else*/ 
  if      (s=="single"     ) emtemp=SINGLE     ;
  else if (s=="ensemble"   ) emtemp=ENSEMBLE   ;
  else if (s=="master"     ) emtemp=MASTER     ;
  else if (s!="master_fast") 
    is.clear(std::ios_base::badbit);

  if (is) em=emtemp;
  return is;
}


evolution::Pars::Pars(parameters::ParameterTable& p, const std::string& mod) 
  : ParsRun(p,mod),
    mcwf::Pars(p,mod),
    evol(p.addTitle("Evolution",mod).addMod("evol",mod,"Evolution mode (single, ensemble, master, master_fast)",SINGLE)),
    negativity(p.addMod("negativity",mod,"Calculates negativity in ensemble & master",false)),
    timeAverage(p.addMod("timeAverage",mod,"Calculates time averages in MCWF trajectory",false)),
    relaxationTime(p.addMod("relaxationTime",mod,"Relaxation time for time averaging",0.))
{}


