#include "ModeFunctionType.h"

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
  case MFT_PLUS : return exp(DCOMP_I*x);
  case MFT_MINUS:                      ;
  }
  return exp(-DCOMP_I*x);
}
