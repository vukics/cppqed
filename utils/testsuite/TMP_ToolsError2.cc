#include "TMP_Tools.h"

void f()
{
  // Should not compile:
  IsEvenAssert<31> error_IsEvenAssert;
}
