#include "BlitzArraySliceIterator.h"

#include "Mode_.h"

#include "impl/LazyDensityOperator.tcc"
#include "StateVector.h"

#include <boost/bind.hpp>

#define BOOST_TEST_MODULE LazyDensityOperator test
#include <boost/test/unit_test.hpp>


using mathutils::fcmp; using mathutils::sqr; using mathutils::sqrAbs; 

using quantumdata::partialTrace;

const double eps=1e-12;

const tmptools::Vector<0> v0=tmptools::Vector<0>();
const tmptools::Vector<1> v1=tmptools::Vector<1>();

const dcomp alpha(2,3), beta(-1,2), gammaa(2,-2);

const mode::StateVector c23(mode::coherent(alpha,78)), c12(mode::coherent(beta,30)), c22(mode::coherent(gammaa,48)), 
  d23(mode::coherent(alpha,78)), d12(mode::coherent(beta,78)), d22(mode::coherent(gammaa,78));

boost::function<double(const mode::LazyDensityOperator&)> photonNumber(static_cast<double(*)(const mode::LazyDensityOperator&)>(mode::photonNumber));

template<int I>
double 
photonNumberRecurse(const quantumdata::LazyDensityOperator<2>& m)
{
  return partialTrace(m,photonNumber,tmptools::Vector<I>(),0.);
}


const mode::DensityOperatorLow
densityOperator(const mode::LazyDensityOperator& m)
{
  size_t dim=m.getDimension();
  mode::DensityOperatorLow res(dim,dim);
  for (int i=0; i<dim; i++) for (int j=0; j<dim; j++) res(i,j)=m(i,j);
  return res;    
}


BOOST_AUTO_TEST_CASE( RANK_TWO )
{
  quantumdata::StateVector<2> psi(c23*c12);

  BOOST_CHECK(!fcmp(photonNumberRecurse<0>(psi),13.,eps));
  BOOST_CHECK(!fcmp(photonNumberRecurse<1>(psi), 5.,eps));

  quantumdata::DensityOperator<2> rho(psi);

  BOOST_CHECK(!fcmp(photonNumberRecurse<0>(rho),13.,eps));
  BOOST_CHECK(!fcmp(photonNumberRecurse<1>(rho), 5.,eps));

  dcomp alpha1(alpha), beta1(beta), alpha2(beta), beta2(gammaa);

  quantumdata::StateVector<2> psii(d23*d12-d12*d22); double norm=psii.renorm();

  BOOST_CHECK(!fcmp((sqrAbs(beta1)+sqrAbs(beta2)-2.*std::real(braket(d23,d12)*braket(d12,d22)*conj(beta1)*beta2))/sqr(norm),
		    photonNumberRecurse<1>(psii),
		    eps)
	      );


  // The most stringent test for the whole of a partial-trace density operator:

  mode::DensityOperator res(-braket(d23,d12)/sqr(norm)*dyad(d22,d12));

  {
    linalg::CMatrix tempView(res.matrixView());
    linalg::calculateTwoTimesRealPartOfSelf(tempView);
  }

  res()+=(d12.dyad()+d22.dyad())/sqr(norm);

  BOOST_CHECK(max(blitzplusplus::sqrAbs(partialTrace(psii,densityOperator,tmptools::Vector<1>(),mode::DensityOperatorLow())-res()))<sqr(eps));

}



BOOST_AUTO_TEST_CASE( RANK_THREE_RECURSIVE )
{
  quantumdata::StateVector<3> psi(c23*c12*c22);

  BOOST_CHECK(!fcmp(partialTrace(psi,photonNumberRecurse<0>,tmptools::Vector<0,1>(),0.),13.,eps));
  BOOST_CHECK(!fcmp(partialTrace(psi,photonNumberRecurse<1>,tmptools::Vector<0,1>(),0.), 5.,eps));
  BOOST_CHECK(!fcmp(partialTrace(psi,photonNumberRecurse<0>,tmptools::Vector<1,2>(),0.), 5.,eps));
  BOOST_CHECK(!fcmp(partialTrace(psi,photonNumberRecurse<1>,tmptools::Vector<1,2>(),0.), 8.,eps));
  BOOST_CHECK(!fcmp(partialTrace(psi,photonNumberRecurse<0>,tmptools::Vector<2,0>(),0.), 8.,eps));
  BOOST_CHECK(!fcmp(partialTrace(psi,photonNumberRecurse<1>,tmptools::Vector<2,0>(),0.),13.,eps));

  // BOOST_CHECK(!fcmp(partialTrace(psi,photonNumberRecurse<0>,v1,0.), 5.,eps));
  // BOOST_CHECK(!fcmp(partialTrace(psi,photonNumberRecurse,tmptools::Vector<2,0>(),0.), 5.,eps));

}



const mode::StateVector e0(mode::fock(2,3)), e1(mode::fock(3,7)), e2(mode::fock(0,1));


int specialChecker(const quantumdata::LazyDensityOperator<3>& m, size_t i1, size_t i2, size_t i3)
{
  BOOST_CHECK(all(m.getDimensions()==TTD_EXTTINY(3)(i1,i2,i3)));
  partialTrace(m,photonNumberRecurse<0>,tmptools::Vector<1,2>(),0.);
  // This is only for checking whether the special implementations are also recursive.
  return 0;
}


BOOST_AUTO_TEST_CASE( THE_SPECIAL_CASE )
{
  quantumdata::StateVector<3> psi(e0*e1*e2);
  partialTrace(psi,bind(specialChecker,_1,3,7,1),tmptools::Vector<0,1,2>(),0); 
  partialTrace(psi,bind(specialChecker,_1,3,1,7),tmptools::Vector<0,2,1>(),0); 

  quantumdata::DensityOperator<3> rho(psi);
  partialTrace(psi,bind(specialChecker,_1,3,7,1),tmptools::Vector<0,1,2>(),0); 
  partialTrace(psi,bind(specialChecker,_1,3,1,7),tmptools::Vector<0,2,1>(),0); 

}
