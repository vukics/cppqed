#include "VectorFromMatrixSliceIterator.tcc"

using namespace blitzplusplus;

template class basi::Iterator<8,tmptools::Range<7,0>,true>;
template class basi::Iterator<8,tmptools::Ordinals<5>,true>;
template class basi::Iterator<8,vfmsi::LeftRight<4,vfmsi::Left>,true>;
