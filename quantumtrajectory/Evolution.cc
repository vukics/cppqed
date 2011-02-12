#include "Evolution.h"

#include<iostream>


std::ostream& operator<<(std::ostream& os, EvolutionMode em)
{
  switch (em) {
    // case EM_CONVERGENCE: return os<<"convergence";
  case EM_SINGLE      : return os<<"single"      ;
  case EM_ENSEMBLE    : return os<<"ensemble"    ;
  case EM_MASTER      : return os<<"master"      ;
  case EM_MASTER_FAST :        os<<"master_fast" ;
  }
  return os;
}


std::istream& operator>>(std::istream& is, EvolutionMode& em) 
{
  EvolutionMode emtemp=EM_MASTER_FAST;
  std::string s;

  is>>s;
  /*if      (s=="convergence") emtemp=EM_CONVERGENCE;
    else*/ 
  if      (s=="single"     ) emtemp=EM_SINGLE     ;
  else if (s=="ensemble"   ) emtemp=EM_ENSEMBLE   ;
  else if (s=="master"     ) emtemp=EM_MASTER     ;
  else if (s!="master_fast") 
    is.clear(std::ios_base::badbit);

  if (is) em=emtemp;
  return is;
}



