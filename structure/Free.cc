#include "Free.h"


namespace structure {


Free::Free(size_t dimension, const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs) 
  : QuantumSystem<1>(dimension), DynamicsBase(realFreqs,complexFreqs) 
{
}


std::ostream& Free::displayMoreParameters(std::ostream& os) const
{
  return DynamicsBase::displayMoreParameters(os<<"# Dimension: "<<getDimension()<<std::endl);
}


} // structure
