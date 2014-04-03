// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED

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


#endif // CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED
