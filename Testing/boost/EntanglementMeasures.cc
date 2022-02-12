// Copyright András Vukics 2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#ifdef  EIGEN3_FOUND

#include "Qbit.h"

#include "EntanglementMeasures.h"

#include <boost/test/unit_test.hpp>

#include <numbers>

template <int i>
constexpr auto tmpVec = tmptools::Vector<i>{};

typedef quantumdata::DensityOperator<2> DO2;
// typedef quantumdata::LazyDensityOperator<3> LDO3;

using quantumdata::negPT;

using namespace qbit;

// example from https://en.wikipedia.org/wiki/Peres%E2%80%93Horodecki_criterion#Example
auto wernerState(double p)
{
  return p*DO2{(state0()*state1()-state1()*state0())/std::numbers::sqrt2} 
    + (0.25*(1-p))*(DO2{state0()*state0()}+DO2{state0()*state1()}+DO2{state1()*state0()}+DO2{state1()*state1()});
}


BOOST_AUTO_TEST_CASE( NEGATIVITY_OF_PARTIAL_TRANSPOSE , * boost::unit_test::tolerance(1e-14) )
{
  {
    DO2 rho{(state0()*state1()+state1()*state0())/std::numbers::sqrt2};
    
    std::cerr<<rho.getArray()<<std::endl<<Eigen::ComplexEigenSolver<Eigen::MatrixX<dcomp>>{
      Eigen::Map<Eigen::MatrixX<dcomp>>{rho.getArray().data(),long(rho.getTotalDimension()),long(rho.getTotalDimension())},false}.eigenvalues()<<std::endl;
      
    auto neg=negPT(rho, tmpVec<0>);
    
    std::cerr<<neg<<std::endl;
    
    BOOST_TEST ( neg == 0.5 ) ;
    
    neg=negPT(rho, tmpVec<1>);
    
    BOOST_TEST ( neg == 0.5 ) ;
  }

  {
    DO2 rho{(state0()+state1())*(state0()+state1())/2.};
    
    std::cerr<<rho.getArray()<<std::endl<<Eigen::ComplexEigenSolver<Eigen::MatrixX<dcomp>>{
      Eigen::Map<Eigen::MatrixX<dcomp>>{rho.getArray().data(),long(rho.getTotalDimension()),long(rho.getTotalDimension())},false}.eigenvalues()<<std::endl;
      
    auto neg=negPT(rho, tmpVec<0>);
    
    std::cerr<<neg<<std::endl;
    
    BOOST_TEST ( neg == 0. ) ;
    
    neg=negPT(rho, tmpVec<1>);
    
    BOOST_TEST ( neg == 0. ) ;
  }

  {
    std::cerr<<wernerState(0.5).getArray()<<std::endl;
    
    BOOST_TEST ( negPT(wernerState(.99),tmpVec<0>) > 0. ) ;

    BOOST_TEST ( negPT(wernerState(.5),tmpVec<0>) > 0. ) ;

    BOOST_TEST ( negPT(wernerState(.34),tmpVec<0>) > 0. ) ;

    BOOST_TEST ( negPT(wernerState(.32),tmpVec<0>) == 0. ) ;

    BOOST_TEST ( negPT(wernerState(.1),tmpVec<0>) == 0. ) ;

  }
  
}


#endif // EIGEN3_FOUND
