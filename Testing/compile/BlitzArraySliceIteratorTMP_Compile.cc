// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzArraySliceIterator.tcc"

using namespace blitzplusplus;

template class basi::Iterator<8,tmptools::Vector<3,1,6>,true>;
template class basi::Iterator<8,tmptools::Range<6,1>,true>;
template class basi::Iterator<8,tmptools::Ordinals<5>,true>;

#include "VectorFromMatrixSliceIterator.tcc"

template class basi::Iterator<8,vfmsi::LeftRight<4,vfmsi::Left>,true>;
template class basi::Iterator<8,vfmsi::LeftRight<4,vfmsi::Right>,true>;
