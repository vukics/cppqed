#include "ExampleMode.h"


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


double photonNumber(const LazyDensityOperator& matrix)
{
  double res=0;
  for (size_t n=1; n<matrix.getDimension(); ++n)
    res+=n*matrix(n);
  return res;
}


double aJumpProba   (const LazyDensityOperator& matrix, double kappa_nPlus1)
{
  return 2.*kappa_nPlus1*photonNumber(matrix);
}

double aDagJumpProba(const LazyDensityOperator& matrix, double kappa_n     )
{
  return 2.*kappa_n*(photonNumber(matrix)+1.);
}


const Tridiagonal aop(size_t dim)
{
  typedef Tridiagonal::Diagonal Diagonal;
  Diagonal diagonal(dim-1);
  return Tridiagonal(Diagonal(),1,Diagonal(),diagonal=sqrt(blitz::tensor::i+1.));
}


const Tridiagonal nop(size_t dim)
{
  Tridiagonal::Diagonal diagonal(dim);
  return Tridiagonal(diagonal=blitz::tensor::i);
}


const Frequencies freqs(double delta, size_t dim)
{
  Tridiagonal::Diagonal diagonal(dim-1); diagonal=dcomp(0,-delta);
  return Frequencies(1,diagonal);
}
