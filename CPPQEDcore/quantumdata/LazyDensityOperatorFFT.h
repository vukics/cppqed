// -*- C++ -*-
#ifndef QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED

#include "LazyDensityOperatorFwd.h"

#include "CMatrix.h"
#include "Exception.h"

#include "FFT.h"

#include "StateVector.h"

#include <boost/mpl/integral_c.hpp>
#include <boost/shared_ptr.hpp>


namespace quantumdata {


void ffTransform(linalg::CVector&, fft::Direction);
void ffTransform(linalg::CMatrix&, fft::Direction);


struct LazyDensityOperatorFFT_NotImplementedException : cpputils::Exception {};


template<typename V, int RANK>
const boost::shared_ptr<const LazyDensityOperator<RANK> > ffTransform(const LazyDensityOperator<RANK>&, fft::Direction);


} // quantumdata


#endif // QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED
