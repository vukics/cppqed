#include "Free.h"


namespace structure {


Free::Free(size_t dimension, const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs) 
  : QuantumSystem<1>(dimension), DynamicsBase(realFreqs,complexFreqs) 
{
}


void Free::displayMoreParameters(std::ostream& os) const
{
  os<<"# Dimension: "<<getDimension()<<std::endl;
  DynamicsBase::displayMoreParameters(os);
}


} // structure
