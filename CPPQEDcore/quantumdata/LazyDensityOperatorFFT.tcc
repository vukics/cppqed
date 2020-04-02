// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_INCLUDED

#include "LazyDensityOperatorFFT.h"

#include "BlitzArraySliceIterator.tcc"
#include "DensityOperator.tcc"
#include "StateVector.tcc"
#include "TMP_Tools.h"

#include <boost/range/algorithm/for_each.hpp>
#include <boost/mpl/for_each.hpp>


namespace quantumdata {


namespace details {


template<int RANK>
class fftWorkerSV
{
public:
  typedef std::shared_ptr<StateVector<RANK> > SV;

  fftWorkerSV(SV psi, fft::Direction dir) : psi_(psi), dir_(dir) {};

  template<int IDX> void operator()(mpl::integral_c<int,IDX>)
  {
    boost::for_each(blitzplusplus::basi::fullRange<tmptools::Vector<IDX> >(psi_->getArray()),
                    boost::bind(static_cast<void(*)(linalg::CVector&, fft::Direction)>(&quantumdata::ffTransform),_1,dir_));
  }

private:
  const SV psi_;
  const fft::Direction dir_;

};


template<int RANK>
class fftWorkerDO
{
public:
  typedef std::shared_ptr<DensityOperator<RANK> > DO;
  
  fftWorkerDO(DO rho, fft::Direction dir) : rho_(rho), dir_(dir) {};

  template<int IDX> void operator()(mpl::integral_c<int,IDX>)
  {
    boost::for_each(blitzplusplus::basi::fullRange<tmptools::Vector<IDX,IDX+RANK> >(rho_->getArray()),
                    boost::bind(static_cast<void(*)(linalg::CMatrix&, fft::Direction)>(&quantumdata::ffTransform),_1,dir_));
  }
  
private:
  const DO rho_;
  const fft::Direction dir_;
  
};


} // details
  
  
template<typename V, int RANK>
const std::shared_ptr<const LazyDensityOperator<RANK> > ffTransform(const LazyDensityOperator<RANK>& matrix, fft::Direction dir)
{
  typedef StateVector    <RANK> SV;
  typedef DensityOperator<RANK> DO;
  
  if      (const auto psi=dynamic_cast<const SV*>(&matrix) ) {
    const std::shared_ptr<SV> res(std::make_shared<SV>(*psi));
    boost::mpl::for_each<V>(details::fftWorkerSV<RANK>(res,dir));
    return res;
  }
  else if (const auto rho=dynamic_cast<const DO*>(&matrix) ) {
    const std::shared_ptr<DO> res(std::make_shared<DO>(*rho));
    boost::mpl::for_each<V>(details::fftWorkerDO<RANK>(res,dir));
    return res;
  }
  else throw LazyDensityOperatorFFT_NotImplementedException();

}


} // quantumdata


#endif // CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_TCC_INCLUDED
