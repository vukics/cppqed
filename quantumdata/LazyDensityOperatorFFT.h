// -*- C++ -*-
#ifndef QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED

#include "LazyDensityOperatorFwd.h"

#include "CMatrix.h"
#include "Exception.h"

#include "FFTFwd.h"

#include "StateVector.h"

#include <boost/mpl/integral_c.hpp>
#include <boost/shared_ptr.hpp>


namespace quantumdata {


void ffTransform(linalg::CVector&, fft::Direction);
void ffTransform(linalg::CMatrix&, fft::Direction);


struct LazyDensityOperatorFFT_NotImplementedException : cpputils::Exception {};


template<typename V, int RANK>
const boost::shared_ptr<const LazyDensityOperator<RANK> > ffTransform(const LazyDensityOperator<RANK>&, fft::Direction);

namespace details {

template<int RANK>
class fftWorkerSV
{
public:
  typedef StateVector<RANK> SV;
  fftWorkerSV(const boost::shared_ptr<SV> psi, fft::Direction dir) : psi_(psi), dir_(dir) {};
  template<int IDX> void operator()(mpl::integral_c<int,IDX> &);
private:
  boost::shared_ptr<SV> psi_;
  fft::Direction dir_;
};

template<int RANK>
class fftWorkerDO
{
public:
  typedef DensityOperator<RANK> DO;
  fftWorkerDO(const boost::shared_ptr<DO> rho, fft::Direction dir) : rho_(rho), dir_(dir) {};
  template<int IDX> void operator()(mpl::integral_c<int,IDX> &);
private:
  boost::shared_ptr<DO> rho_;
  fft::Direction dir_;
};

} // details

} // quantumdata


#endif // QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED
