// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Mode.tcc"

#include "TridiagonalHamiltonian.h"

#include <boost/assign/list_inserter.hpp>


using std::cout; using std::endl; using std::string;
using namespace boost;
using namespace cppqedutils;


namespace mode {

#define DEFINE_make_by_redirect(AUX) const Ptr make(const BOOST_PP_CAT(Pars,AUX) & p, QM_Picture qmp) {return make<Averaged>(p,qmp);}

DEFINE_make_by_redirect()
DEFINE_make_by_redirect(Lossy)
DEFINE_make_by_redirect(Pumped)
DEFINE_make_by_redirect(PumpedLossy)

#undef DEFINE_make_by_redirect


const Tridiagonal aop(size_t);


////////
//
// Exact
//
////////


Exact::Exact(dcomp zI, double omegaKerr, size_t dim)
  : FreeExact(dim), zI_(zI), omegaKerr_(omegaKerr)
{
}


void Exact::updateU(Time t) const
{
  getDiagonal()=exp(-zI_*(double(t)*blitz::tensor::i));
}


//////////////
//
// Hamiltonian
//
//////////////

#define KERRDIAGONAL(d,o) d=(i-1.)*i; d*=o/DCOMP_I;

const Tridiagonal::Diagonal mainDiagonal(dcomp z, double omegaKerr, size_t dim)
{
  using blitz::tensor::i;
  Tridiagonal::Diagonal res1(dim), res2(dim);
  res1=i; res1*=z;
  KERRDIAGONAL(res2,omegaKerr)
  return Tridiagonal::Diagonal(res1+res2);
}


const Tridiagonal pumping(dcomp eta, size_t dim)
{
  return tridiagMinusHC(eta*aop(dim).dagger());
}


namespace {


Tridiagonal hOverI(dcomp z, dcomp eta, double omegaKerr, size_t dim)
// Here we use runtime dispatching because runtime dispatching will anyway happen at the latest in some maker function.
{
  bool isPumped=isNonZero(eta);
  if (isNonZero(z)) {
    Tridiagonal res(mainDiagonal(-z,omegaKerr,dim));
    if (isPumped) return res+pumping(eta,dim);
    else return res;
  }
  else if(isPumped) return pumping(eta,dim);
  else return Tridiagonal();
}


}


#define FILLKERR using blitz::tensor::i; Tridiagonal::Diagonal temp(dim); KERRDIAGONAL(temp,omegaKerrAlter) getH_OverIs().push_back(Tridiagonal(temp));


template<> template<>
Hamiltonian<true >::Hamiltonian(dcomp zSch, dcomp zI, dcomp eta, double omegaKerr, double omegaKerrAlter, size_t dim)
  : Base(furnishWithFreqs(hOverI(zSch,eta,omegaKerr,dim),mainDiagonal(zI,omegaKerr,dim))), Exact(zI,omegaKerr,dim), zSch_(zSch), eta_(eta), dim_(dim)
{
  FILLKERR
}

template<> template<>
Hamiltonian<false>::Hamiltonian(dcomp zSch, dcomp eta, double omegaKerr, double omegaKerrAlter, size_t dim)
  : Base(hOverI(zSch,eta,omegaKerr,dim)), zSch_(zSch), eta_(eta), dim_(dim)
{
  FILLKERR
}

#undef KERRDIAGONAL
#undef FILLKERR

//////////////
//
// Liouvillean
//
//////////////

void details::aJump   (StateVectorLow& psi, double kappa) // kappa is kappa*(nTh+1) when the temperature is finite
{
  double fact=sqrt(2.*kappa);
  int ubound=psi.ubound(0);
  for (int n=0; n<ubound; ++n)
    psi(n)=fact*sqrt(n+1)*psi(n+1);
  psi(ubound)=0;
}

void details::aDagJump(StateVectorLow& psi, double kappa) // kappa is kappa* nTh    when the temperature is finite
{
  double fact=sqrt(2.*kappa);
  for (int n=psi.ubound(0); n>0; --n)
    psi(n)=fact*sqrt(n)*psi(n-1);
  psi(0)=0;
}

void details::aSuperoperator(const DensityOperatorLow& rho, DensityOperatorLow& drhodt, double kappa)
{
  for (int m=0; m<rho.ubound(0); ++m)
    for (int n=0; n<rho.ubound(1); ++n)
      drhodt(m,n)+=2*kappa*sqrt((m+1)*(n+1))*rho(m+1,n+1);
}

void details::aDagSuperoperator(const DensityOperatorLow& rho, DensityOperatorLow& drhodt, double kappa)
{
  for (int m=1; m<=rho.ubound(0); ++m)
    for (int n=1; n<=rho.ubound(1); ++n)
      drhodt(m,n)+=2*kappa*sqrt(m*n)*rho(m-1,n-1);
}


///////////
//
// Averaged
//
///////////


namespace {

typedef cppqedutils::KeyPrinter::KeyLabels KeyLabels;

const KeyLabels Assemble(const KeyLabels& first, const KeyLabels& middle, const KeyLabels& last=KeyLabels())
{
  KeyLabels res(first); boost::assign::push_back(res).range(middle).range(last);
  return res;
}

}


Averaged::Averaged(const KeyLabels& follow, const KeyLabels& precede)
  : Base(keyTitle,
         Assemble(precede,KeyLabels{"<number operator>","VAR(number operator)","real(<ladder operator>)","imag(\")"},follow))
{
}


const Averages Averaged::average_v(NoTime, const LazyDensityOperator& matrix) const
{
  auto averages(initializedAverages());

  for (int n=1; n<int(matrix.getDimension()); n++) {

    double diag=matrix(n);
    averages(0)+=  n*diag;
    averages(1)+=n*n*diag;

    dcomp offdiag(sqrt(n)*matrix(n)(n-1));
    averages(2)+=real(offdiag);
    averages(3)+=imag(offdiag);

  }

  return averages;

}


void Averaged::process_v(Averages& averages) const
{
  averages(1)-=sqr(averages(0));
}



AveragedQuadratures::AveragedQuadratures(const KeyLabels& follow, const KeyLabels& precede)
  : Averaged(Assemble(KeyLabels{"VAR(X)","VAR(Y)","COV(X,Y)"},follow),precede)
{
}


const Averages AveragedQuadratures::average_v(NoTime t, const LazyDensityOperator& matrix) const
{
  auto averages(initializedAverages());

  averages(blitz::Range(0,3))=Averaged::average_v(t,matrix)(blitz::Range(0,3));

  for (int n=2; n<int(matrix.getDimension()); n++) {

    dcomp  offdiag(sqrt(n*(n-1.))*matrix(n)(n-2));
    averages(4)+=real(offdiag);
    averages(5)+=imag(offdiag);

  }

  return averages;

}


void AveragedQuadratures::process_v(Averages& averages) const
{
  {
    Averages ranged(averages(blitz::Range(0,3)));
    Averaged::process_v(ranged);
  }

  double 
    Re_aSqrExpVal=averages(4),
    Re_aExpVal=averages(2),
    Im_aExpVal=averages(3);
    
  double
    a=averages(0)+.5+Re_aSqrExpVal-2*sqr(Re_aExpVal), // <X^2>-<X>^2
    b=averages(0)+.5-Re_aSqrExpVal-2*sqr(Im_aExpVal), // <Y^2>-<Y>^2
    c=averages(5)-2*Re_aExpVal*Im_aExpVal; // <(XY+YX)/2>-<X><Y>

  averages(4)=a;
  averages(5)=b;
  averages(6)=c;
  
  /*
  Eigenvalues of the real matrix
    ( a c )
    ( c b )
  */

  /*
    double
    temp1=(a+b)/2,
    temp2=sqrt(sqr(a-b)+4*sqr(c))/2;

    averages(4)=temp1+temp2;
    averages(5)=temp1-temp2;
    
    averages(6)=atan((-a+b)/(2*c)+temp2/c);
  */
}


//////////
//
// Helpers
//
//////////


const Tridiagonal aop(size_t dim)
{
  typedef Tridiagonal::Diagonal Diagonal;
  Diagonal diagonal(dim-1);
  return Tridiagonal(Diagonal(),1,Diagonal(),diagonal=sqrt(blitz::tensor::i+1.));
}



// This returns a Tridiagonal furnished with frequencies, when mode is derived from mode::Exact

const Tridiagonal aop(Ptr mode)
{
  size_t dim=mode->getDimension();
  Tridiagonal res(aop(dim));
  if (const auto exact=dynamic_cast<const mode::Exact*>(mode.get())) res.furnishWithFreqs(mainDiagonal(exact->get_zI(),exact->get_omegaKerr(),dim));
  return res;
}



const Tridiagonal nop(Ptr mode)
{
  return Tridiagonal(mainDiagonal(1.,0.,mode->getDimension()));
}



double photonNumber(const StateVectorLow& psi)
{
  using blitz::tensor::i;
  return sum(i*blitzplusplus::sqrAbs(psi(i)));
}


double photonNumber(const LazyDensityOperator& matrix)
{
  double res=0;
  for (size_t n=1; n<matrix.getDimension(); ++n)
    res+=n*matrix(n);
  return res;
}


StateVector coherent(dcomp alpha, size_t dim)
{
  StateVector res(dim,false);
  double norm(exp(-sqr(abs(alpha))/2.));

  for (size_t n=0; n<dim; ++n) res(n)=norm*coherentElement(n,alpha);

  return res;

}


StateVector fock(size_t n, size_t dim, double phase)
{
  if (n>=dim) throw std::overflow_error("Fock state "+std::to_string(n)+" higher than dim "+std::to_string(dim));
  StateVector res(dim);
  res(n)=exp(DCOMP_I*phase);
  return res;
}


StateVector init(const Pars& p)
{
  return p.minitFock ? fock(p.minitFock,p.cutoff) : coherent(p.minit,p.cutoff);
}


} // mode




ModeBase::ModeBase(size_t dim, const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs, const string& keyTitle)
    : Free(dim,realFreqs,complexFreqs)
{
  getParsStream()<<keyTitle<<endl;
}




PumpedLossyModeIP_NoExact::PumpedLossyModeIP_NoExact(const mode::ParsPumpedLossy& p)
  : ModeBase(p.cutoff),
    quantumoperator::TridiagonalHamiltonian<1,true>(furnishWithFreqs(mode::pumping(p.eta,p.cutoff),mode::mainDiagonal(dcomp(p.kappa,-p.delta),p.omegaKerr,p.cutoff))),
    structure::ElementLiouvillean<1,1,true>(mode::keyTitle,"excitation loss"),
    structure::ElementAveraged<1,true>(mode::keyTitle,{"<number operator>","real(<ladder operator>)","imag(\")"}),
    z_(p.kappa,-p.delta)
{
  getParsStream()<<"Interaction picture, not derived from Exact\n";
}



double PumpedLossyModeIP_NoExact::rate(OneTime, const LazyDensityOperator& m) const
{
  return mode::photonNumber(m);
}


void PumpedLossyModeIP_NoExact::doActWithJ(OneTime t, StateVectorLow& psi) const
{
  dcomp fact=sqrt(2.*real(z_))*exp(-z_*double(t));
  int ubound=psi.ubound(0);
  for (int n=0; n<ubound; ++n)
    psi(n)=fact*sqrt(n+1)*psi(n+1);
  psi(ubound)=0;
}



const mode::Averages PumpedLossyModeIP_NoExact::average_v(OneTime t, const LazyDensityOperator& matrix) const
{
  auto averages(initializedAverages());

  for (int n=1; n<int(matrix.getDimension()); n++) {

    averages(0)+=n*exp(-2*real(z_)*n*t)*matrix(n);

    dcomp offdiag(sqrt(n)*matrix(n)(n-1)*exp(-(z_+2*(n-1)*real(z_))*double(t)));
    averages(1)+=real(offdiag);
    averages(2)+=imag(offdiag);

  }

  return averages;

}
