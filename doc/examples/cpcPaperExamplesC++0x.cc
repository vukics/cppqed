#include <blitz/array.h>

//#include "Types.h"
/*
#include "LazyDensityOperator.h"
#include "StateVector.h"

#include "BlitzArraySliceIterator.h"
#include "Range.h"

#include <valarray>

*/
template <int RANK> using StateVector=blitz::Array<std::complex<double>,RANK> quantumdata::Types<RANK>::StateVectorLow;

/*
template <int RANK>
struct DensityOperator : quantumdata::Types<RANK>::DensityOperatorLow
{
  typedef typename quantumdata::Types<RANK>::DensityOperatorLow Base;
  
  DensityOperator(Base b) : Base(b) {}

  DensityOperator(size_t d) : Base(d) {}

};


using tmptools::Vector;

#define BASI_Range blitzplusplus::basi::fullRange

using blitzplusplus::basi::begin;
using cpputils::for_each;

using quantumdata::LazyDensityOperator;

typedef dcomp complex;
*/

// ****************************************
// ****************************************
// ****************************************


void actWithA(StateVector<5>&);


void actOnExtended(StateVector<11>& psi)
{
  for_each(BASI_Range(psi, 
                      Vector<3,6,1,9,7>()),
           actWithA);
}


void TwoModesParticle2D_Hamiltonian(const StateVector<4>&, StateVector<4>&);


void Hamiltonian(const StateVector<8>& psi, StateVector<8>& dpsidt)
{
  for_each(BASI_Range(psi, Vector<4,0,7,1>()), 
           begin(dpsidt, Vector<4,0,7,1>()),
           TwoModesParticle2D_Hamiltonian);
}


void composeWithA(DensityOperator<5>& rho)
{
  for_each(BASI_Range(rho, Vector<0,1,2,3,4>()),
           actWithA);
}


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
  typedef LazyDensityOperator<2>::Idx Idx; 
  // The index type is a blitz::TinyVector<ptrdiff_t,2>
  const LazyDensityOperator<2>::Dimensions dim(m.getDimensions());

  complex res;
  for (int i=0; i<dim[0]-1; ++i) for (int j=1; j<dim[1]; ++j)
    res+=sqrt((i+1)*j)*matrix(Idx(i,j),Idx(i+1,j-1));
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
