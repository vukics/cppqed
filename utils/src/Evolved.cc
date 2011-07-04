#include "Evolved.h"

namespace evolved {

  
details::EvolvedCommon::EvolvedCommon(double dtInit, double epsRel, double epsAbs)
  : t_(0), dtTry_(dtInit), dtDid_(0), epsRel_(epsRel), epsAbs_(epsAbs) 
{}




} // evolved
