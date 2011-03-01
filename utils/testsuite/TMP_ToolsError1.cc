#include "TMP_Tools.h"

void f()
{
  // Should not compile:
  typedef tmptools::OrdinalMF<-1>::type ordinalError;
}
