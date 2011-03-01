#include "TMP_Tools.h"

void f()
{
  // Should not compile:
  tmptools::pair_c<23,23>();
  /*
  P_23_42::SanityCheck<24,45>();
  P_23_42::SanityCheck<21,40>();
  */
}
