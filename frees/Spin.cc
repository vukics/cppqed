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


const Tridiagonal splus(size_t twos)
{
  typedef Tridiagonal::Diagonal Diagonal;
  Diagonal diagonal(twos);
  using blitz::tensor::i;
  return Tridiagonal(Diagonal(),1,diagonal=sqrt((twos-i)*(i+1)));
}


const Tridiagonal sz(size_t twos)
{
  Tridiagonal::Diagonal diagonal(twos+1);
  return Tridiagonal(diagonal=blitz::tensor::i-twos/2.);
}



const Frequencies freqs(const SpinBase* spin)
{
  if (dynamic_cast<const structure::FreeExact*>(spin)) return mode::freqs(spin->get_z(),spin->getDimension());
  else return Frequencies(spin->getDimension());
}


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : twos(p.addMod<size_t>("twos",mod,"2*s, the size of the spin (dimension twos+1)",1)),
    omega(p.addMod("omega",mod,"Spin precession frequency",1.)),
    gamma(p.addMod("gamma",mod,"Spin decay rate"          ,1.))
{}



} // spin


SpinBase::SpinBase(size_t twos, double omega, double gamma) 
  : Free(twos+1,tuple_list_of("omega",omega,1)("gamma",gamma,1)), structure::ElementAveraged<1>("Spin",list_of("<sz>")("<sz^2>")("real(<s^+>)")("imag(\")")), omega_(omega), gamma_(gamma)
{
  getParsStream()<<"# Spin "<<twos/2.<<endl;
}



const SpinBase::Averages SpinBase::average(const LazyDensityOperator& matrix) const
{
  double s=(getDimension()-1.)/2.;

  Averages averages(4);

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
