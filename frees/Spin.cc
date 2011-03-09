#include "Spin.h"

#include "Mode.h"

#include "LazyDensityOperator.h"

#include "Pars.h"

#include<boost/assign/list_of.hpp>
#include<boost/assign/std/list.hpp>

using namespace std;
using namespace boost;
using namespace assign;
using namespace mathutils;

namespace spin {

size_t decideDimension(size_t twos, size_t dim)
{
  return dim && dim<twos+1 ? dim : twos+1;
}


const Tridiagonal splus(size_t twos, size_t dim)
{
  typedef Tridiagonal::Diagonal Diagonal;
  Diagonal diagonal(decideDimension(twos,dim)-1);
  using blitz::tensor::i;
  return Tridiagonal(Diagonal(),1,diagonal=sqrt((twos-i)*(i+1)));
}


const Tridiagonal sz(size_t twos, size_t dim)
{
  Tridiagonal::Diagonal diagonal(decideDimension(twos,dim));
  return Tridiagonal(diagonal=blitz::tensor::i-twos/2.);
}



const Frequencies freqs(const SpinBase* spin)
{
  if (dynamic_cast<const structure::FreeExact*>(spin)) return mode::freqs(spin->get_z(),spin->getDimension());
  else return Frequencies(spin->getDimension());
}


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : twos(p.addTitle("Spin",mod).addMod<size_t>("twos",mod,"2*s, the size of the spin (dimension: twos+1 or spinDim)",1)),
    dim (p.addMod<size_t>("spinDim",mod,"the dimension of the truncated spin Hilbert space",0)),
    omega(p.addMod("omega",mod,"Spin precession frequency",1.)),
    gamma(p.addMod("gamma",mod,"Spin decay rate"          ,1.))
{}



} // spin


SpinBase::SpinBase(size_t twos, double omega, double gamma, size_t dim) 
  : Free(spin::decideDimension(twos,dim),tuple_list_of("omega",omega,1)("gamma",gamma,1)),
    structure::ElementAveraged<1>("Spin",list_of("<sz>")("<sz^2>")("real(<s^+>)")("imag(\")")("|Psi(dim-1)|^2")), omega_(omega), gamma_(gamma)
{
  getParsStream()<<"# Spin "<<twos/2.<<endl;
}



const SpinBase::Averages SpinBase::average(const LazyDensityOperator& matrix) const
{
  double s=(getDimension()-1.)/2.;

  Averages averages(5);

  averages=0;

  for (int n=0; n<matrix.getDimension(); n++) {

    double diag=matrix(n);
    averages(0)+=   (n-s)*diag;
    averages(1)+=mathutils::sqr(n-s)*diag;

    if (n<2*s) {
      dcomp temp(sqrt((2*s-n)*(n+1))*matrix(n,n+1));
      averages(2)+=real(temp);
      averages(3)+=imag(temp);
    }

  }

  averages(4)=matrix(matrix.getDimension()-1);

  return averages;

}


void SpinBase::process(Averages& averages) const
{
  averages(1)-=sqr(averages(0));
}


void Spin::updateU(double dtdid) const
{
  Factors& factors(getFactors());
  for (int i=0; i<factors.size(); i++)
    factors(i)=exp(-dtdid*i*get_z());
}
