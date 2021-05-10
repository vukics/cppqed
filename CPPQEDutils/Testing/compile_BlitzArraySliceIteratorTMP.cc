// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzArray.h"
#include "SliceIterator.h"

template class cppqedutils::SliceIterator<CArray,8,tmptools::Vector<3,1,6>>;
template class cppqedutils::SliceIterator<CArray,8,tmptools::Range<1,7>>;
template class cppqedutils::SliceIterator<CArray,8,tmptools::CopyToVector<tmptools::Ordinals<5>>>;

/*
#include "VectorFromMatrixSliceIterator.h"

namespace blitzplusplus { namespace vfmsi {

template class Iterator<8,Left >;
template class Iterator<8,Right>;

} }
*/
