#include "Mode.h"

#include "ParsMode.h"

#include "StateVector.h"

#include<boost/bind.hpp>

#include<boost/preprocessor/control/if.hpp>
#include<boost/preprocessor/comparison/equal.hpp>

using namespace std;
using namespace boost;
using namespace assign;
using namespace mathutils;


namespace mode {


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


namespace details {


const Tridiagonal pumping(const dcomp& eta, size_t dim)
{
  return tridiagMinusHC(eta*aop(dim).dagger());
}


const Tridiagonal hOverI(const dcomp& z, const dcomp& eta, size_t dim)
// Here we use runtime dispatching because runtime dispatching will
// anyway happen at the latest in some maker function.
{
  bool isPumped=isNonZero(eta);
  if (isNonZero(z)) {
    Tridiagonal res(-z*nop(dim));
    if (isPumped) return res+pumping(eta,dim);
    else return res;
  }
  else if(isPumped) return pumping(eta,dim);
  else return Tridiagonal();
}


} // details


template<>
Hamiltonian<true >::Hamiltonian(const dcomp& zSch, const dcomp& zI, const dcomp& eta, size_t dim, boost::mpl::true_)
  : Base(details::hOverI(zSch,eta,dim),freqs(zI,dim)), Exact(zI,dim), zSch_(zSch), eta_(eta), dim_(dim)
{}

template<>
Hamiltonian<false>::Hamiltonian(const dcomp& zSch, const dcomp& eta, size_t dim, boost::mpl::false_)
  : Base(details::hOverI(zSch,eta,dim)), zSch_(zSch), eta_(eta), dim_(dim)
{}




//////////////
//
// Liouvillean
//
//////////////


namespace details {


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


} // details


Liouvillean<true >::Liouvillean(const ParsLossy& p) 
  : Base(JumpStrategies(bind(details::aJump   ,_1,p.kappa*(p.nTh+1)),
			bind(details::aDagJump,_1,p.kappa* p.nTh  )),
	 JumpProbabilityStrategies(bind(details::aJumpProba   ,_1,p.kappa*(p.nTh+1)),
				   bind(details::aDagJumpProba,_1,p.kappa* p.nTh  ))),
    kappa_(p.kappa), nTh_(p.nTh)
{
}

template<>
void Liouvillean<false,false>::doActWithJ(StateVectorLow& psi) const
{
  details::aJump(psi,kappa_);
}

template<>
void Liouvillean<false,true>::doActWithJ(StateVectorLow& psi) const
{
  details::aJump(psi,kappa_);
}


template<>
double Liouvillean<false,false>::probability(const LazyDensityOperator& matrix) const
{
  return details::aJumpProba(matrix,kappa_);
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

typedef structure::ElementAveragedCommon::KeyLabels KeyLabels;

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
  : Base("Mode",
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

    double sqrtn=sqrt(double(n));
    dcomp offdiag(matrix(n,n-1));
    averages(2)+=sqrtn*real(offdiag);
    averages(3)+=sqrtn*imag(offdiag);

  }

  return averages;

}


void Averaged::process(Averages& averages) const
{
  averages(1)-=sqr(averages(0));
}



AveragedQuadratures::AveragedQuadratures(const KeyLabels& follow, const KeyLabels& precede)
  : Averaged(Assemble(KeyLabels(),list_of("<X^2>-<X>^2")("<Y^2>-<Y>^2")("<(XY+YX)/2>-<X><Y>"),follow),precede)
    // Last parameter is necessary, otherwise ambiguity with Averaged copy constructor
{
}


const AveragedQuadratures::Averages AveragedQuadratures::average(const LazyDensityOperator& matrix) const
{
  Averages averages(7);

  averages=0;

  averages(blitz::Range(0,3))=Averaged::average(matrix);

  for (int n=2; n<int(matrix.getDimension()); n++) {

    double sqrtn=sqrt(double(n));
    dcomp  offdiag(matrix(n,n-2));

    averages(4)+=sqrtn*sqrt(n-1.)*real(offdiag);
    averages(5)+=sqrtn*sqrt(n-1.)*imag(offdiag);

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


const Tridiagonal aop(const ModeBase* mode)
{
  return aop(mode->getDimension());
}


const Tridiagonal nop(size_t dim)
{
  Tridiagonal::Diagonal diagonal(dim);
  return Tridiagonal(diagonal=blitz::tensor::i);
}


const Tridiagonal nop(const ModeBase* mode)
{
  return nop(mode->getDimension());
}


const Tridiagonal xop(size_t dim)
{
  return tridiagPlusHC(aop(dim))/sqrt(2.);
}


const Tridiagonal xop(const ModeBase* mode)
{
  return xop(mode->getDimension());
}


const Frequencies freqs(const dcomp& zI, size_t dim)
{
  Tridiagonal::Diagonal diagonal(dim-1); diagonal=zI;
  return Frequencies(1,diagonal);
}


const Frequencies freqs(const ModeBase* mode)
{
  if (const mode::Exact* exact=dynamic_cast<const mode::Exact*>(mode)) return freqs(exact->get_zI(),mode->getDimension());
  else return Frequencies(mode->getDimension());
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




ModeBase::ModeBase(size_t dim, const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs)
    : Free(dim,realFreqs,complexFreqs)
{
  getParsStream()<<"# Mode\n";
}


