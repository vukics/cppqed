// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionMethod.h"

#include<iostream>


std::ostream& evolution::operator<<(std::ostream& os, Method em)
{
  switch (em) {
    // case CONVERGENCE: return os<<"convergence";
  case SINGLE   : return os<<"single"  ;
  case ENSEMBLE : return os<<"ensemble";
  case MASTER   :        os<<"master"  ;
  }
  return os;
}


std::istream& evolution::operator>>(std::istream& is, Method& em) 
{
  Method emtemp=MASTER;
  std::string s;

  is>>s;
  /*if      (s=="convergence") emtemp=CONVERGENCE;
    else*/ 
  if      (s=="single"  ) emtemp=SINGLE     ;
  else if (s=="ensemble") emtemp=ENSEMBLE   ;
  else if (s!="master") 
    is.clear(std::ios_base::badbit);

  if (is) em=emtemp;
  return is;
}


