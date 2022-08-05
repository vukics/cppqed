// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Spin.h"

#include "MathExtensions.h"

#include "Mode_.h"

#include "LazyDensityOperator.h"
#include "Tridiagonal.tcc"

#include "Pars.h"

#include<boost/assign/std/list.hpp>

using namespace std;
using namespace boost;
using namespace assign;
using namespace cppqedutils;

namespace spin {

namespace {

typedef Tridiagonal::Diagonal Diagonal;

}


size_t decideDimension(size_t twoS, size_t dim)
{
  return (dim && dim<twoS+1) ? dim : twoS+1;
}


Diagonal mainDiagonal(size_t dim, size_t twoS, dcomp z)
{
  Diagonal diagonal(dim);
  diagonal=blitz::tensor::i-twoS/2.;
  return Diagonal(z*diagonal);
}


Tridiagonal splus(size_t dim, size_t twoS, dcomp z, bool isExact)
{
  Diagonal diagonal(dim-1);
  using blitz::tensor::i;
  // i=m+s (m is the magnetic quantum number)
  Tridiagonal res(Diagonal(),1,diagonal=sqrt((twoS-i)*(i+1.)));
  // writing "1." is tremendously important here, for converting the whole thing to doubles: otherwise, the sqrt is apparently performed within integers by blitz!!! 
  if (isExact) res.furnishWithFreqs(mainDiagonal(dim,twoS,z));
  return res;
}


Tridiagonal sn(Ptr spin)
{
  return sn(spin->getDimension(),spin->getTwoS(),spin->get_z(),bool(dynamic_pointer_cast<const structure::FreeExact<false>>(spin)),spin->getTheta(),spin->getPhi());
} 


Tridiagonal sz(size_t dim, size_t twoS)
{
  Diagonal diagonal(dim);
  return Tridiagonal(diagonal=blitz::tensor::i-twoS/2.);
}

Tridiagonal sz(Ptr spin) {return sz(spin->getDimension(),spin->getTwoS());}


Tridiagonal splus(Ptr spin) {return splus(spin->getDimension(),spin->getTwoS(),spin->get_z(),bool(dynamic_pointer_cast<const structure::FreeExact<false>>(spin)));}



Pars::Pars(parameters::Table& p, const std::string& mod)
  : twoS(p.addTitle("Spin",mod).add<size_t>("twoS",mod,"2*s, the size of the spin (dimension: twoS+1 or spinDim)",1)),
    dim (p.add<size_t>("spinDim",mod,"the dimension of the truncated spin Hilbert space",0)),
    theta(p.add("theta",mod,"Spin orientation inclination",0.)),
    phi  (p.add("phi"  ,mod,"Spin orientation azimuth    ",0.)),
    omega(p.add("omega",mod,"Spin precession frequency",1.)),
    gamma(p.add("gamma",mod,"Spin decay rate"          ,1.))
{}



} // spin


SpinBase::SpinBase(size_t twoS, double theta, double phi, double omega, double gamma, size_t dim) 
  : Free(spin::decideDimension(twoS,dim),{RF{"omega",omega,1},RF{"gamma",gamma,1}}),
    structure::ElementAveraged<1>("Spin",{"<sz>","<sz^2>","real(<s^+>)","imag(\")","|Psi(dim-1)|^2"}), twoS_(twoS), theta_(theta), phi_(phi), omega_(omega), gamma_(gamma), s_(twoS/2.)
{
  getParsStream()<<"Spin "<<s_<<endl;
  getParsStream()<<"theta="<<theta_<<" phi="<<phi_<<endl;
}



const spin::Averages SpinBase::average_v(structure::NoTime, const LazyDensityOperator& matrix) const
{
  auto averages(initializedAverages());

  for (int n=0; n<matrix.getDimension(); n++) {

    double diag=matrix(n);
    averages(0)+=(n-s_)*diag;
    averages(1)+=sqr(n-s_)*diag;

    if (n<matrix.getDimension()-1) {
      dcomp temp(sqrt((twoS_-n)*(n+1))*matrix(n)(n+1));
      averages(2)+=real(temp);
      averages(3)+=imag(temp);
    }

  }

  averages(4)=matrix(matrix.getDimension()-1);

  return averages;

}


void SpinBase::process_v(spin::Averages& averages) const
{
  averages(1)-=sqr(averages(0));
}


void Spin::updateU(structure::OneTime dtDid) const
{
  Diagonal& factors(getDiagonal());
  for (int i=0; i<factors.size(); i++)
    factors(i)=exp(-dtDid*i*get_z());
}


void spin::Liouvillian::doActWithJ(structure::NoTime, structure::freesystem::StateVectorLow& psi) const
{
  double fact=sqrt(2.*gamma_);
  int ubound=psi.ubound(0);
  for (int n=0; n<ubound; ++n)
    psi(n)=fact*sqrt((n+1)*(twoS_-n))*psi(n+1);
  psi(ubound)=0;
}


double spin::Liouvillian::rate(structure::NoTime, const structure::freesystem::LazyDensityOperator& matrix) const
{
  double res=0;
  for (size_t n=1; n<matrix.getDimension(); ++n)
    res+=2.*gamma_*n*(twoS_-n+1)*matrix(n);
  return res;
}
