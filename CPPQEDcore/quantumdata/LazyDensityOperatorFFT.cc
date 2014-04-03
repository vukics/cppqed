// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "LazyDensityOperatorFFT.h"

#include "BlitzArrayTraits.h"
#include "FFT.tcc"
#include "VectorFromMatrixSliceIterator.tcc"

#include <boost/range/algorithm/for_each.hpp>


using namespace fft; using namespace linalg;


namespace {

  
void ffTransformCV(CVector& psi, Direction dir)
{
  struct Helper
  {
    static void _(CVector& psi, int i1, int i2, double norm)
    {
      dcomp temp(psi(i1));
      psi(i1)=norm*psi(i2);
      psi(i2)=norm*temp;
    }
  };

  int size=psi.size();

  if (size<2) return;

  transform(psi,dir);

  int halfnumber=psi.size()>>1;

  // NEEDS_WORK express the following with blitz
  for (int j=0; j<halfnumber; j++ ) Helper::_(psi,j,j+halfnumber,pow(size,-.5));
  for (int j=1; j<size      ; j+=2) psi(j)*=-1;

}


}


void quantumdata::ffTransform(linalg::CVector& psi, fft::Direction dir)
{
  ffTransformCV(psi,dir);
}



void quantumdata::ffTransform(linalg::CMatrix& rho, fft::Direction dir)
{
  using namespace blitzplusplus::vfmsi;

  for_each(fullRange<Left >(rho),bind(&ffTransformCV,_1,        dir ));
  for_each(fullRange<Right>(rho),bind(&ffTransformCV,_1,reverse(dir)));

}


