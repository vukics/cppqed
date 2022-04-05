// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED

#include "DensityOperator.h"
#include "StateVector.h"

#include "CMatrix.h"
#include "FFT.h"
#include "SliceIterator.h"

#include <boost/range/algorithm/for_each.hpp>

#include <memory>


namespace mpl=boost::mpl;


namespace quantumdata {


void ffTransform(linalg::CVector&, fft::Direction);
void ffTransform(linalg::CMatrix&, fft::Direction);



template<typename V, int RANK>
const std::shared_ptr<const LazyDensityOperator<RANK> > ffTransform(const LazyDensityOperator<RANK>& matrix, fft::Direction dir, V v=V{})
{
  typedef StateVector    <RANK> SV;
  typedef DensityOperator<RANK> DO;
  
  if      (const auto psi=dynamic_cast<const SV*>(&matrix) ) {
    const std::shared_ptr<SV> res(std::make_shared<SV>(*psi));
    hana::for_each(v,[&](auto a) {
      boost::for_each(cppqedutils::sliceiterator::fullRange<tmptools::Vector<a.value> >(psi->getArray()),
                      [=](linalg::CVector& psiS){ffTransform(psiS,dir);});
    });
    return res;
  }
  else if (const auto rho=dynamic_cast<const DO*>(&matrix) ) {
    const std::shared_ptr<DO> res(std::make_shared<DO>(*rho));
    hana::for_each(v,[&](auto a) {
      constexpr auto IDX=a.value;
      boost::for_each(cppqedutils::sliceiterator::fullRange<tmptools::Vector<IDX,IDX+RANK> >(rho->getArray()),
                      [=](linalg::CMatrix& rhoS){ffTransform(rhoS,dir);});
    });
    return res;
  }
  else throw std::runtime_error("LazyDensityOperator FFT not implemented");

}


} // quantumdata


#endif // CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORFFT_H_INCLUDED
