// Copyright András Vukics 2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#ifdef  EIGEN3_FOUND

#include "Qbit.h"

#include "EntanglementMeasures.h"

#include <boost/test/unit_test.hpp>

#include <numbers>


template <int i>
constexpr auto tmpVec = tmptools::Vector<i>{};

using DO1=quantumdata::DensityOperator<1>;
using DO2=quantumdata::DensityOperator<2>;

// typedef quantumdata::LazyDensityOperator<3> LDO3;

using quantumdata::negPT; using quantumdata::mutualInformation;

using namespace qbit;

// example from https://en.wikipedia.org/wiki/Peres%E2%80%93Horodecki_criterion#Example
DO2 wernerState(double p)
{
  return p*DO2{(state0()*state1()-state1()*state0())/std::numbers::sqrt2} 
    + (0.25*(1-p))*(DO2{state0()*state0()}+DO2{state0()*state1()}+DO2{state1()*state0()}+DO2{state1()*state1()});
}


const DO2 rhoBell{(state0()*state1()+state1()*state0())/std::numbers::sqrt2};


BOOST_AUTO_TEST_CASE( NEGATIVITY_OF_PARTIAL_TRANSPOSE , * boost::unit_test::tolerance(1e-14) )
{
  {
    std::cerr<<rhoBell.getArray()<<std::endl<<Eigen::ComplexEigenSolver<Eigen::MatrixX<dcomp>>{
      Eigen::Map<Eigen::MatrixX<dcomp>>{const_cast<dcomp*>(rhoBell.getArray().data()),long(rhoBell.getTotalDimension()),long(rhoBell.getTotalDimension())},false
    }.eigenvalues()<<std::endl;
      
    auto neg=negPT(rhoBell, tmpVec<0>);
    
    std::cerr<<neg<<std::endl;
    
    BOOST_TEST ( neg == 0.5 ) ;
    
    neg=negPT(rhoBell, tmpVec<1>);
    
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


BOOST_AUTO_TEST_CASE( ENTROPY_AND_MUTUAL_INFORMATION , * boost::unit_test::tolerance(1e-15) )
{
  static_assert( mpl::equal<tmptools::NegatedView<10,tmptools::Vector<0,5,1,3> >, tmptools::Vector<2,4,6,7,8,9> >::value ) ;
 
  /* reduce<tmptools::Vector<0>>(wernerState(0)); reduce<tmptools::NegatedVector<2,tmptools::Vector<0>>>(wernerState(0));*/
  
  DO1 rhoReduced0{reduce<tmptools::Vector<0>>(rhoBell)} , rhoReduced1{reduce<tmptools::Vector<1>>(rhoBell)} , rhoMixed{{2},true};
  rhoMixed(0)(0)=.5; rhoMixed(1)(1)=.5;
  
  BOOST_TEST ( frobeniusNorm( rhoReduced0-rhoMixed ) == 0 ) ;

  BOOST_TEST ( frobeniusNorm( rhoReduced1-rhoMixed ) == 0 ) ;
  
  BOOST_TEST ( mutualInformation(rhoBell,tmpVec<0>) == -2.*std::log(.5) ) ;
  
  BOOST_TEST ( mutualInformation(wernerState(0),tmpVec<0>) == 0. ) ;
  
  BOOST_TEST ( mutualInformation(wernerState(1),tmpVec<0>) == -2.*std::log(.5) ) ;
  
  BOOST_TEST ( mutualInformation(wernerState(.5),tmpVec<0>) == 0.3127515147113673 ) ;
  
  /*
  std::cerr.precision(16);
  std::cerr<<mutualInformation(rhoBell,tmpVec<0>)<<std::endl;
  
  std::cerr<<mutualInformation(wernerState(0),tmpVec<0>)<<std::endl
    <<mutualInformation(wernerState(1),tmpVec<0>)<<std::endl
    <<mutualInformation(wernerState(.5),tmpVec<0>)<<std::endl;
  */
}

#endif // EIGEN3_FOUND
