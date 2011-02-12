#include "Frequencies.h"

namespace quantumoperator {


template<>
Frequencies<1>::Frequencies(size_t size, mpl::int_<1>)
  : freqs_(blitzplusplus::TOA_ShallowCopy(),Diagonal(size),Diagonal(),Diagonal())
{
  freqs_(0)=0.;
}


template<>
Frequencies<1>::Frequencies(size_t K, const Diagonal& minus, mpl::int_<1>)
  : freqs_(blitzplusplus::TOA_DeepCopy(),Diagonal(minus.size()+K),minus,Diagonal(-minus))
{
  freqs_(0)=0.;
}


/*
const Frequencies<1> zeroFreq()
{
  return Frequencies<1>(0);
}
*/

} // quantumoperator
