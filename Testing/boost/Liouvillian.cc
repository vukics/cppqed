#include "Qbit.h"

#include <boost/test/unit_test.hpp>

using namespace cppqedutils; using namespace quantumdata;

BOOST_AUTO_TEST_CASE( QBIT_RATE_AND_SUPEROPERATOR , * boost::unit_test::tolerance(1e-15) )
{
  using namespace qbit;

  StateVector<1> psi({2}); psi(0)={2.,-3.}; psi(1)={-1.,2.}; renorm(psi);
  double gamma=1.;

  BOOST_TEST( structure::rateFromJump<1>(0.,psi,sigmaJump(gamma)) == sigmaRate(gamma)(psi) ) ;
  BOOST_TEST( structure::rateFromJump<1>(0.,psi,sigmaPlusJump(gamma)) == sigmaPlusRate(gamma)(psi) );

  const DensityOperator<1> rho(psi);
  const auto rowIterationOffsets{calculateSlicesOffsets<retainedAxes<0>>(rho.extents)};

  {
    DensityOperator<1> drhodtDirect{{2},zeroInit<2>}, drhodtIndirect{{2},zeroInit<2>};

    sigmaSuperoperator(gamma)(rho,drhodtDirect.mutableView());

    superoperatorFromJump<1>(0.,rho,drhodtIndirect,sigmaJump(gamma),rowIterationOffsets);

    std::cerr<<json(drhodtDirect)<<std::endl<<json(drhodtIndirect)<<std::endl;

    BOOST_TEST( frobeniusNorm(drhodtDirect-drhodtIndirect) == 0 );
  }

  {
    DensityOperator<1> drhodtDirect{{2},zeroInit<2>}, drhodtIndirect{{2},zeroInit<2>};

    sigmaPlusSuperoperator(gamma)(rho,drhodtDirect.mutableView());

    superoperatorFromJump<1>(0.,rho,drhodtIndirect,sigmaPlusJump(gamma),rowIterationOffsets);

    std::cerr<<json(drhodtDirect)<<std::endl<<json(drhodtIndirect)<<std::endl;

    BOOST_TEST( frobeniusNorm(drhodtDirect-drhodtIndirect) == 0 );
  }

  // StateVector<2> psi2(psi*psi);

}

