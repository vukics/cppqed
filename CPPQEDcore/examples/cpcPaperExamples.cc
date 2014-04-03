// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "VectorFromMatrixSliceIterator.h"

#include <boost/range/algorithm/for_each.hpp>

using boost::for_each; using tmptools::Vector;

template <int RANK> using StateVector    =CArray<  RANK>;
template <int RANK> using DensityOperator=CArray<2*RANK>;


void actWithA(StateVector<5>&);


void actOnExtended(StateVector<11>& psi)
{
  using namespace blitzplusplus::basi;
  for_each(fullRange<Vector<3,6,1,9,7> >(psi),actWithA);
}


void composeWithA(DensityOperator<5>& rho)
{
  using namespace blitzplusplus::vfmsi;
  for_each(fullRange<Left>(rho),actWithA);
}

/*
const complex
calculateASqr(const LazyDensityOperator<1>& m)
{
  complex res;
  for (int i=2; i<m.getTotalDimension(); ++i)
    res+=sqrt(i*(i-2))*m(i,i-2);
  return res;
}


const complex
calculateADaggerB(const LazyDensityOperator<2>& m)
{
  // The index type is a blitz::TinyVector<ptrdiff_t,2>
  const LazyDensityOperator<2>::Dimensions dim(m.getDimensions());

  complex res;
  for (int i=0; i<dim[0]-1; ++i) for (int j=1; j<dim[1]; ++j)
    res+=sqrt((i+1)*j)*matrix(i,j)(i+1,j-1);
  return res;
}


template <int RANK, typename F, typename V, typename T>
const T
partialTrace(const LazyDensityOperator<RANK>& matrix, F f, V v, T t);


const DensityOperator<1>
copyDensityOperator(const LazyDensityOperator<1>& m)
{
  size_t dim=m.getDimension();
  DensityOperator<1> res(dim);
  for (int i=0; i<dim; i++) for (int j=0; j<dim; j++)
    res(i,j)=m(i,j);
  return res;
}


template<int RANK, int SUBSYSTEM>
// for the subsystem defined by the index SUBSYSTEM
const DensityOperator<1>
partialTraceOfUnarySubsystem(const LazyDensityOperator<RANK>& m)
{
  return partialTrace(m,
                      copyDensityOperator,
                      Vector<SUBSYSTEM>(),
                      DensityOperator<1>());
}


template<int RANK>
// RANK>3
const complex
calculateADaggerB_atPositions3and1(const LazyDensityOperator<RANK>& m)
{
  return partialTrace(m,
                      calculateADaggerB,
                      Vector<3,1>(),
                      complex());
}


typedef std::valarray<double> Averages;

const Averages
calculateModeAverages(const LazyDensityOperator<1>& m)
{
  Averages averages(0.,4);

  for (int n=1; n<m.getDimension(); n++) {
    double diag=m(n);
    averages[0]+=  n*diag;
    averages[1]+=n*n*diag;

    double sqrtn=sqrt(double(n)); complex offdiag(m(n,n-1));
    averages[2]+=sqrtn*real(offdiag);
    averages[3]+=sqrtn*imag(offdiag);
  }

  return averages;

}


template<int RANK, int MODE_POSITION>
// for the subsystem defined by the index MODE_POSITION
const Averages
calculateEmbeddedModeAverages(const LazyDensityOperator<RANK>& m)
{
  return partialTrace(m,
                      calculateModeAverages,
                      Vector<MODE_POSITION>(),
                      Averages());
}


// ****************************************
// ****************************************
// ****************************************

// Explicit instantiations:

quantumdata::StateVector<4> psi(43);

template const complex calculateADaggerB_atPositions3and1(const LazyDensityOperator<4>& m);

complex c=calculateADaggerB_atPositions3and1(psi);

template const Averages calculateEmbeddedModeAverages<4,3>(const LazyDensityOperator<4>& m);
*/