// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ModeFunction.h"

#include<iostream>


using namespace std;

ostream& operator<<(ostream& os, ModeFunctionType mf)
{
  switch (mf) {
  case MFT_SIN  : return os<<"Sin"  ;
  case MFT_COS  : return os<<"Cos"  ;
  case MFT_PLUS : return os<<"Plus" ;
  case MFT_MINUS:        os<<"Minus";
  }
  return os;
}

istream& operator>>(istream& is, ModeFunctionType& mf) 
{
  ModeFunctionType mftemp=MFT_MINUS;
  string s;

  is>>s;
       if (s=="Sin"  ) mftemp=MFT_SIN ;
  else if (s=="Cos"  ) mftemp=MFT_COS ;
  else if (s=="Plus" ) mftemp=MFT_PLUS;
  else if (s!="Minus") 
    is.clear(ios_base::badbit);

  if (is) mf=mftemp;
  return is;
}



const dcomp modeFunction(ModeFunctionType mf, double x)
{
  switch (mf) {
  case MFT_SIN  : return sin(x)        ;
  case MFT_COS  : return cos(x)        ;
  case MFT_PLUS : return exp(1i*x);
  case MFT_MINUS:                      ;
  }
  return exp(-1i*x);
}


/*
std::ostream& operator<<(std::ostream& os, const ModeFunction& mf)
{
  os<<"mode-function type is "<<mf.get<0>()<<", wave number: "<<mf.get<1>();
  return os;
}


*/
