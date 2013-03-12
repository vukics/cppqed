// -*- C++ -*-
#ifndef QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_H_INCLUDED

#include "LazyDensityOperatorFFT.h"

#include "impl/BlitzArraySliceIterator.tcc"
#include "impl/DensityOperator.tcc"
#include "impl/StateVector.tcc"
#include "TMP_Tools.h"

#include <boost/range/algorithm/for_each.hpp>
#include <boost/make_shared.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/vector.hpp>

namespace quantumdata {


namespace details {

template<int RANK>
template<int IDX>
void fftWorkerSV<RANK>::operator()(mpl::integral_c<int,IDX> &)
{
  using tmptools::Vector;
  using blitzplusplus::basi::fullRange;
  void (* const ffTransformWorker)(linalg::CVector&, fft::Direction) = &quantumdata::ffTransform;
  boost::for_each(fullRange<Vector<IDX> >((*psi_)()),boost::bind(ffTransformWorker,_1,dir_));
}

template<int RANK>
template<int IDX>
void fftWorkerDO<RANK>::operator()(mpl::integral_c<int,IDX> &)
{
  using tmptools::Vector;
  using blitzplusplus::basi::fullRange;
  void (* const ffTransformWorker)(linalg::CMatrix&, fft::Direction) = &quantumdata::ffTransform;
  boost::for_each(fullRange<Vector<IDX,IDX+RANK> >((*rho_)()),boost::bind(ffTransformWorker,_1,dir_));
}

} // details
  
  
template<typename V, int RANK>
const boost::shared_ptr<const LazyDensityOperator<RANK> > ffTransform(const LazyDensityOperator<RANK>& matrix, fft::Direction dir)
{
  typedef StateVector    <RANK> SV;
  typedef DensityOperator<RANK> DO;
  
  
  if      (const SV*const psi=dynamic_cast<const SV*>(&matrix) ) {
    const boost::shared_ptr<SV> res(boost::make_shared<SV>(*psi));
    boost::mpl::for_each<V>(details::fftWorkerSV<RANK>(res,dir));
    return res;
  }
  else if (const DO*const rho=dynamic_cast<const DO*>(&matrix) ) {
    const boost::shared_ptr<DO> res(boost::make_shared<DO>(*rho));
    boost::mpl::for_each<V>(details::fftWorkerDO<RANK>(res,dir));
    return res;
  }
  else throw LazyDensityOperatorFFT_NotImplementedException();

}


} // quantumdata


#endif // QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_H_INCLUDED
