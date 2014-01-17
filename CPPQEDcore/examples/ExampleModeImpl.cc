#include "ExampleMode.h"
#include "Tridiagonal.tcc"

void aJump   (StateVectorLow& psi, double kappa_nPlus1)
{
  double fact=sqrt(2.*kappa_nPlus1);
  int ubound=psi.ubound(0);
  for (int n=0; n<ubound; ++n)
    psi(n)=fact*sqrt(n+1)*psi(n+1);
  psi(ubound)=0;
}


void aDagJump(StateVectorLow& psi, double kappa_n     )
{
  double fact=sqrt(2.*kappa_n);
  for (int n=psi.ubound(0); n>0; --n)
    psi(n)=fact*sqrt(n)*psi(n-1);
  psi(0)=0;
}


namespace {

double photonNumber(const LazyDensityOperator& matrix)
{
  double res=0;
  for (size_t n=1; n<matrix.getDimension(); ++n)
    res+=n*matrix(n);
  return res;
}

}


double aJumpRate   (const LazyDensityOperator& matrix, double kappa_nPlus1)
{
  return 2.*kappa_nPlus1*photonNumber(matrix);
}


double aDagJumpRate(const LazyDensityOperator& matrix, double kappa_n     )
{
  return 2.*kappa_n*(photonNumber(matrix)+matrix.trace());
}


const Tridiagonal aop(size_t dim)
{
  typedef Tridiagonal::Diagonal Diagonal;
  Diagonal diagonal(dim-1);
  return Tridiagonal(Diagonal(),1,Diagonal(),diagonal=sqrt(blitz::tensor::i+1.));
}


const Tridiagonal::Diagonal mainDiagonal(const dcomp& z, size_t dim)
{
  Tridiagonal::Diagonal res(dim);
  res=blitz::tensor::i;
  return res*=z;
}


const Tridiagonal nop(size_t dim)
{
  return Tridiagonal(mainDiagonal(1.,dim));
}


const Tridiagonal aop(const hierarchical::ModeBase& mode)
{
  using namespace hierarchical;
  size_t dim=mode.getDimension();
  if (const auto modeIP=dynamic_cast<const PumpedLossyModeIP*>(&mode))
    return furnishWithFreqs(aop(dim),mainDiagonal(modeIP->get_z(),dim));
  else
    return aop(dim);
}
