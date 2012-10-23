#include "FreeExact.h"

namespace structure {

using namespace std;

void FreeExact::actWithU_v(double dtDid, StateVectorLow& psi) const
{
  if (dtDid!=dtDid_) {updateU(dtDid_=dtDid);}
  psi*=factors_;
}

} // Structure
