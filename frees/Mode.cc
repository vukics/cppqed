#include "Mode_.h"

#include "ParsMode.h"

#include "impl/StateVector.tcc"

#include "impl/TridiagonalHamiltonian.tcc"

#include <boost/assign.hpp>
#include <boost/bind.hpp>


using namespace std;
using namespace boost;
using namespace assign;
using namespace mathutils;


namespace mode {

const Tridiagonal aop(size_t);


////////
//
// Exact
//
////////


Exact::Exact(const dcomp& zI, size_t dim)
  : FreeExact(dim), zI_(zI)
{
}


void Exact::updateU(double dtdid) const
{
  getFactors()=exp(-zI_*(dtdid*blitz::tensor::i));
}


//////////////
//
// Hamiltonian
//
//////////////


const Tridiagonal::Diagonal mainDiagonal(const dcomp& z, size_t dim)
{
  Tridiagonal::Diagonal res(dim);
  res=blitz::tensor::i;
  return res*=z;
}


const Tridiagonal pumping(const dcomp& eta, size_t dim)
{
  return tridiagMinusHC(eta*aop(dim).dagger());
}


namespace {


const Tridiagonal hOverI(const dcomp& z, const dcomp& eta, size_t dim)
// Here we use runtime dispatching because runtime dispatching will
// anyway happen at the latest in some maker function.
{
  bool isPumped=isNonZero(eta);
  if (isNonZero(z)) {
    Tridiagonal res(mainDiagonal(-z,dim));
    if (isPumped) return res+pumping(eta,dim);
    else return res;
  }
  else if(isPumped) return pumping(eta,dim);
  else return Tridiagonal();
}


}


template<>
Hamiltonian<true >::Hamiltonian(const dcomp& zSch, const dcomp& zI, const dcomp& eta, size_t dim, boost::mpl::true_)
  : Base(furnishWithFreqs(hOverI(zSch,eta,dim),mainDiagonal(zI,dim))), Exact(zI,dim), zSch_(zSch), eta_(eta), dim_(dim)
{}

template<>
Hamiltonian<false>::Hamiltonian(const dcomp& zSch, const dcomp& eta, size_t dim, boost::mpl::false_)
  : Base(hOverI(zSch,eta,dim)), zSch_(zSch), eta_(eta), dim_(dim)
{}




//////////////
//
// Liouvillean
//
//////////////

namespace {


void aJump   (StateVectorLow& psi, double kappa)
{
  double fact=sqrt(2.*kappa);
  int ubound=psi.ubound(0);
  for (int n=0; n<ubound; ++n)
    psi(n)=fact*sqrt(n+1)*psi(n+1);
  psi(ubound)=0;
}


void aDagJump(StateVectorLow& psi, double kappa)
{
  double fact=sqrt(2.*kappa);
  for (int n=psi.ubound(0); n>0; --n)
    psi(n)=fact*sqrt(n)*psi(n-1);
  psi(0)=0;
}


double aJumpProba   (const LazyDensityOperator& matrix, double kappa)
{
  return 2.*kappa*photonNumber(matrix);
}


double aDagJumpProba(const LazyDensityOperator& matrix, double kappa)
{
  return 2.*kappa*(photonNumber(matrix)+1.);
}


}


Liouvillean<true >::Liouvillean(double kappa, double nTh, const std::string& kT)
  : Base(JumpStrategies(bind(aJump   ,_1,kappa*(nTh+1)),
			bind(aDagJump,_1,kappa* nTh  )),
	 JumpProbabilityStrategies(bind(aJumpProba   ,_1,kappa*(nTh+1)),
				   bind(aDagJumpProba,_1,kappa* nTh  )),
	 kT,list_of("excitation loss")("excitation absorption")),
    kappa_(kappa), nTh_(nTh)
{
}

template<>
void Liouvillean<false,false>::doActWithJ(StateVectorLow& psi) const
{
  aJump(psi,kappa_);
}

template<>
void Liouvillean<false,true>::doActWithJ(StateVectorLow& psi) const
{
  aJump(psi,kappa_);
}


template<>
double Liouvillean<false,false>::probability(const LazyDensityOperator& matrix) const
{
  return aJumpProba(matrix,kappa_);
}


template<>
double Liouvillean<false,true>::probability(const LazyDensityOperator&) const
{
  return -1;
}


///////////
//
// Averaged
//
///////////


namespace {

typedef cpputils::KeyPrinter::KeyLabels KeyLabels;

const KeyLabels Assemble(const KeyLabels& first,
			 const KeyLabels& middle,
			 const KeyLabels& last=KeyLabels()
			 )
{
  KeyLabels res(first);
  boost::copy(last,boost::copy(middle,std::back_inserter(res)));
  return res;
}

}


Averaged::Averaged(const KeyLabels& follow, const KeyLabels& precede)
  : Base(keyTitle,
	 Assemble(precede,list_of("<number operator>")("VAR(number operator)")("real(<ladder operator>)")("imag(\")"),follow))
{
}


const Averaged::Averages Averaged::average(const LazyDensityOperator& matrix) const
{
  Averages averages(4);

  averages=0;

  for (int n=1; n<int(matrix.getDimension()); n++) {

    double diag=matrix(n);
    averages(0)+=  n*diag;
    averages(1)+=n*n*diag;

    dcomp offdiag(sqrt(n)*matrix(n,n-1));
    averages(2)+=real(offdiag);
    averages(3)+=imag(offdiag);

  }

  return averages;

}


void Averaged::process(Averages& averages) const
{
  averages(1)-=sqr(averages(0));
}



AveragedQuadratures::AveragedQuadratures(const KeyLabels& follow, const KeyLabels& precede)
  : Averaged(Assemble(KeyLabels(),list_of("VAR(X)")("VAR(Y)")("COV(X,Y)"),follow),precede)
    // Last parameter is necessary, otherwise ambiguity with Averaged copy constructor
{
}


const AveragedQuadratures::Averages AveragedQuadratures::average(const LazyDensityOperator& matrix) const
{
  Averages averages(7);

  averages=0;

  averages(blitz::Range(0,3))=Averaged::average(matrix);

  for (int n=2; n<int(matrix.getDimension()); n++) {

    dcomp  offdiag(sqrt(n*(n-1.))*matrix(n,n-2));
    averages(4)+=real(offdiag);
    averages(5)+=imag(offdiag);

  }

  return averages;

}


void AveragedQuadratures::process(Averages& averages) const
{
  {
    Averages ranged(averages(blitz::Range(0,3)));
    Averaged::process(ranged);
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

const Tridiagonal aop(SmartPtr mode)
{
  size_t dim=mode->getDimension();
  Tridiagonal res(aop(dim));
  if (const mode::Exact* exact=dynamic_cast<const mode::Exact*>(mode.get())) res.furnishWithFreqs(mainDiagonal(exact->get_zI(),dim));
  return res;
}



const Tridiagonal nop(SmartPtr mode)
{
  return Tridiagonal(mainDiagonal(1.,mode->getDimension()));
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


const StateVector coherent(const dcomp& alpha, size_t dim)
{
  bool caught=false;

  StateVector res(dim,false);
  double norm(exp(-sqr(abs(alpha))/2.));

  for (size_t n=0; n<dim; ++n) {
    try {
      dcomp temp(norm*pow(alpha,n)/sqrt(fact(n)));
      res()(n)=temp;
    } catch (FactOverflow) {
      if (!caught) cout<<"# LossyMode beware: Coherent called with too high photon numbers!"<<endl;
      caught=true;
    }
  }

  return res;

}


const StateVector fock(size_t n, size_t dim, double phase) throw(PrepError)
{
  if (n>=dim) throw PrepError();
  StateVector res(dim);
  res()(n)=exp(DCOMP_I*phase);
  return res;
}


const StateVector init(const Pars& p)
{
  return p.minitFock ? fock(p.minitFock,p.cutoff) : coherent(p.minit,p.cutoff);
}


} // mode




ModeBase::ModeBase(size_t dim, const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs, const string& keyTitle)
    : Free(dim,realFreqs,complexFreqs)
{
  getParsStream()<<"# "<<keyTitle<<endl;
}





PumpedLossyModeIP_NoExact::PumpedLossyModeIP_NoExact(const mode::ParsPumpedLossy& p)
  : ModeBase(p.cutoff),
    structure::TridiagonalHamiltonian<1,true>(furnishWithFreqs(mode::pumping(p.eta,p.cutoff),
							       mode::mainDiagonal(dcomp(p.kappa,-p.delta),p.cutoff))),
    structure::ElementLiouvillean<1,1,true>(mode::keyTitle,"excitation loss"),
    structure::ElementAveraged<1,true>(mode::keyTitle,list_of("<number operator>")("real(<ladder operator>)")("imag(\")")),
    z_(p.kappa,-p.delta)
{
  getParsStream()<<"# Interaction picture, not derived from Exact\n";
}



double PumpedLossyModeIP_NoExact::probability(double, const LazyDensityOperator& m) const
{
  return mode::photonNumber(m);
}


void PumpedLossyModeIP_NoExact::doActWithJ(double t, StateVectorLow& psi) const
{
  dcomp fact=sqrt(2.*real(z_))*exp(-z_*t);
  int ubound=psi.ubound(0);
  for (int n=0; n<ubound; ++n)
    psi(n)=fact*sqrt(n+1)*psi(n+1);
  psi(ubound)=0;
}



const PumpedLossyModeIP_NoExact::Averages PumpedLossyModeIP_NoExact::average(double t, const LazyDensityOperator& matrix) const
{
  Averages averages(3);

  averages=0;

  for (int n=1; n<int(matrix.getDimension()); n++) {

    averages(0)+=n*exp(-2*real(z_)*n*t)*matrix(n);

    dcomp offdiag(sqrt(n)*matrix(n,n-1)*exp(-(z_+2*(n-1)*real(z_))*t));
    averages(1)+=real(offdiag);
    averages(2)+=imag(offdiag);

  }

  return averages;

}
