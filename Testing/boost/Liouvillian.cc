#include "Qbit.h"

#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( QBIT_RATE_AND_SUPEROPERATOR , * boost::unit_test::tolerance(1e-15) )
{
  using namespace qbit;

  quantumdata::StateVector<1> psi({2}); psi(0)={2.,-3.}; psi(1)={-1.,2.};
  double gamma=1.;

  BOOST_TEST( structure::rateFromJump<1>(0.,psi,sigmaJump(gamma)) == sigmaRate(gamma)(psi) ) ;
  BOOST_TEST( structure::rateFromJump<1>(0.,psi,sigmaPlusJump(gamma)) == sigmaPlusRate(gamma)(psi) );

  quantumdata::DensityOperator<1> rho(psi); rho=0.;

  // quantumdata::StateVector<2> psi2(psi*psi);

}

