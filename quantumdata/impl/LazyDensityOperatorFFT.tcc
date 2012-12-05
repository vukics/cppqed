// -*- C++ -*-
#ifndef QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_H_INCLUDED

#include "LazyDensityOperatorFFT.h"

#include "impl/DensityOperator.tcc"
#include "impl/StateVector.tcc"

#include <boost/make_shared.hpp>

namespace quantumdata {


namespace details {
  
  
} // details
  
  
template<int RANK>
const boost::shared_ptr<const LazyDensityOperator<RANK> > ffTransform(const LazyDensityOperator<RANK>& matrix, fft::Direction dir)
{
  typedef StateVector    <RANK> SV;
  typedef DensityOperator<RANK> DO;
  
  if      (const SV*const psi=dynamic_cast<const SV*>(&matrix) ) {
    const boost::shared_ptr<SV> res(boost::make_shared<SV>(*psi));
    linalg::CVector view(res->vectorView());
    ffTransform(view,dir);
    return res;
  }
  else if (const DO*const rho=dynamic_cast<const DO*>(&matrix) ) {
    const boost::shared_ptr<DO> res(boost::make_shared<DO>(*rho));
    linalg::CMatrix view(res->matrixView());
    ffTransform(view,dir);
    return res;
  }
  else throw LazyDensityOperatorFFT_NotImplementedException();

}


} // quantumdata


#endif // QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_H_INCLUDED