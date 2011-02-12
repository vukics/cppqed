/*

  let x, k, x0, k0 be 2D vectors and A a real, orthogonal, positive definite matrix

  Then if 

  F(x)=exp(-(x-x0)*A*(x-x0)+i*k0*x)

  then its Fourier transform is

  F(k) \propto exp(-(k-k0)*Ainv*(k-k0)/4-i*x0*k)  

  Where the proportionality factor is to be established.

*/

#include "Particle.h"

#include "Types.h"

#include "BlitzArrayTraits.h"
#include "FFT.h"
// #include "MathExtensions.h"
#include "VectorFromMatrixSliceIterator.h"

#include <flens/flens.h>


using namespace particle;
using namespace flens;
using namespace std;
using namespace blitzplusplus;
using namespace vfmsi;
using namespace fft;

typedef structure::Types<2>::StateVectorLow SVL;

typedef GeMatrix<FullStorage<double, ColMajor> >  GEMatrix;
typedef DenseVector<Array<double> >        DEVector;

int main()
{

  Spatial space(7);
  size_t dim=space.getDimension();

  SVL psi(shape(dim,dim)), psiK(shape(dim,dim));

  DEVector x0(2), k0(2);
  x0=.4, -.5; 
  k0=-20, 10;

  GEMatrix A(2,2);
  A=
    9.8, 1.6, 
    1.6, 8.1;  

  GEMatrix Ainv(A);

  {
    DenseVector<Array<int> > p(2);
    trf(Ainv,p);
    tri(Ainv,p);
  }

  cerr<<A<<endl<<Ainv<<endl;


  for (int i=0; i<dim; i++)
    for (int j=0; j<dim; j++) {

      DEVector x(2), k(2);
      x=space.x(i), space.x(j);
      k=space.k(i), space.k(j);
      double 
	xAx   =DEVector(x-x0)*DEVector(A   *DEVector(x-x0)), 
	kAinvk=DEVector(k-k0)*DEVector(Ainv*DEVector(k-k0));

      psi (i,j)=exp(-xAx      +DCOMP_I*(k0*x));
      psiK(i,j)=exp(-kAinvk/4.-DCOMP_I*(x0*k));

    }

  psiK*=exp(DCOMP_I*(x0*k0))*mathutils::PI/sqrt(A(1,1)*A(2,2)-A(1,2)*A(2,1));

  SVL psiFFT(psi.copy()), psiKFFT(psiK.copy());

  for_each(begin<1,vfmsi::Left ,false>(psiFFT),end<1,vfmsi::Left ,false>(psiFFT),bind(&ffTransform,_1,FFTDIR_XK));
  for_each(begin<1,vfmsi::Right,false>(psiFFT),end<1,vfmsi::Right,false>(psiFFT),bind(&ffTransform,_1,FFTDIR_XK));

  for_each(begin<1,vfmsi::Left ,false>(psiKFFT),end<1,vfmsi::Left ,false>(psiKFFT),bind(&ffTransform,_1,FFTDIR_KX));
  for_each(begin<1,vfmsi::Right,false>(psiKFFT),end<1,vfmsi::Right,false>(psiKFFT),bind(&ffTransform,_1,FFTDIR_KX));

  /*
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) 

      cout<<i<<' '<<j<<' '<<real(psi(i,j))<<' '<<imag(psi(i,j))<<' '<<real(psiK(i,j))<<' '<<imag(psiK(i,j))<<endl;

    cout<<endl;
  }
  */

  cerr<<
    max(abs(SVL(psi-psi(50,60)/psiKFFT(50,60)*psiKFFT)))<<' '<<
    max(abs(SVL(psiK-psiK(40,70)/psiFFT(40,70)*psiFFT)))<<endl<<
    psi(50,60)/psiKFFT(50,60)<<' '<<psiK(40,70)/psiFFT(40,70)<<endl;

  return 0;

}
