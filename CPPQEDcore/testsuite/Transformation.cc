// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#define BOOST_TEST_MODULE Transformation
#include <boost/test/unit_test.hpp>

#include "ComplexArrayExtensions.h"
#include "Transformation.h"

#include "Randomized.h"

#include <boost/function.hpp>

#include <iostream>


using namespace quantumdata;
using namespace blitzplusplus;


typedef CArray<1> Array1;
typedef CArray<2> Array2;
typedef CArray<4> Array4;
typedef CArray<6> Array6;

typedef transformation::Traits<Array2> Array2Traits;
typedef transformation::Traits<Array4> Array4Traits;

typedef transformation::Identity<1> Identity1;
typedef transformation::Traits<Identity1> Identity1Traits;

typedef transformation::Identity<2> Identity2;



template<int RANK>
struct Dummy {};


namespace quantumdata {

namespace transformation {


template<int RANK>
struct Traits<Dummy<RANK> >
{
  static const int N_RANK=RANK;
};


} // transformation

} // quantumdata



BOOST_AUTO_TEST_CASE(CompositeAuxiliaryVector)
{
  typedef boost::mpl::first<transformation::namehider::Algorithm<boost::fusion::result_of::make_vector<Dummy<2>,transformation::Identity<3>,Dummy<1>,Dummy<3>,transformation::Identity<4>,Dummy<3> >::type>::type>::type AuxiliaryVector;

  BOOST_STATIC_ASSERT((tmptools::numerical_equal<AuxiliaryVector,boost::mpl::vector_c<int,0,5,6,13> >::value));

}


BOOST_AUTO_TEST_CASE(LinearTransformation)
{
  Array1 x(3), x_(3);
  x=1,2,3;

  Array2 T(3,3);
  T=
    1,0,0,
    0,1,dcomp(0.5,0.5),
    0,dcomp(0.5,-0.5),1;
  Array2Traits::transform(T, x, x_);
  BOOST_CHECK_EQUAL(x_(0), T(0,0)*x(0) + T(0,1)*x(1) + T(0,2)*x(2));
  BOOST_CHECK_EQUAL(x_(1), T(1,0)*x(0) + T(1,1)*x(1) + T(1,2)*x(2));
  BOOST_CHECK_EQUAL(x_(2), T(2,0)*x(0) + T(2,1)*x(1) + T(2,2)*x(2));
}

BOOST_AUTO_TEST_CASE(IdentityTransformation)
{
  Array1 x(3);
  x=1,2,3;
  Identity1 I;
  Array1 y(3);
  Identity1Traits::transform(I, x, y);

  BOOST_CHECK(all(y==x));
}


BOOST_AUTO_TEST_CASE(CompositeTransformation)
{
  Array1 x(2), x_(2), y(3), y_(3);

  x=1,2; y=3,4,-1;

  Array2 T1(2,2), T2(3,3);
  T1=
    1,dcomp(0,1),
    dcomp(0,-1),1;
  T2=
     2, 1,-1,
     1, 3,-2,
    -1,-2, 5;

  Array2 xy(2,3), xy_1(2,3), xy_2(2,3);
  typedef transformation::Compose<Array2,Array2> Composite;
  Composite::type c(Composite::compose(T1,T2));
  typedef transformation::Traits<Composite::type> CTraits;

  Array2Traits::transform(T1,x,x_);
  Array2Traits::transform(T2,y,y_);

  xy_1 = blitzplusplus::doDirect<blitzplusplus::dodirect::multiplication>(x_,y_);

  xy   = blitzplusplus::doDirect<blitzplusplus::dodirect::multiplication>(x ,y );

  CTraits::transform(c,xy,xy_2);

  BOOST_CHECK(blitz::all(xy_1==xy_2));

  // std::cout<<max(sqrAbs(xy_1-xy_2));

}


void trafoFunction(const Array2& in, Array2& out) {out=in; swap(out(1,2),out(2,1));}


BOOST_AUTO_TEST_CASE(MoreCompositeTransformation)
{
  using randomized::fillWithRandom; using namespace blitzplusplus;

  Array1 
    psi1(6), psi1v(psi1.shape()),
    psi4(4), psi4v(psi4.shape()); // x_(2), y(3), y_(3);

  Array2 
    psi23(5,3), psi23v(psi23.shape()),
    psi56(3,2), psi56v(psi56.shape());

  Array2 t1(6,6), t4(4,4);

  Array4 t56(3,2,3,2);

  fillWithRandom(t56,fillWithRandom(t4,fillWithRandom(t1,fillWithRandom(psi56,fillWithRandom(psi4,fillWithRandom(psi23,fillWithRandom(psi1)))))));

  t1 +=hermitianConjugate(t1 );
  t4 +=hermitianConjugate(t4 );
  t56+=hermitianConjugate(t56);

  Array6 psiAll(6,5,3,4,3,2), psiAllv(psiAll.shape()), psiAllvv(psiAll.shape());

  ////////////////////////////////////
  // Composite trafo with one Identity
  ////////////////////////////////////

  Array2Traits::transform(t1,psi1,psi1v);
  Array2Traits::transform(t4,psi4,psi4v);

  Array4Traits::transform(t56,psi56,psi56v);

  {
    using namespace blitz::tensor;
    psiAll =psi1 (i)*psi23(j,k)*psi4 (l)*psi56 (m,n);
    psiAllv=psi1v(i)*psi23(j,k)*psi4v(l)*psi56v(m,n);
  }

  {
    using namespace boost::fusion;
    transformation::Composite<result_of::make_vector<Array2,Identity2,Array2,Array4>::type>
      (make_vector(t1,Identity2(),t4,t56))
      .transform(psiAll,psiAllvv);
  }

  BOOST_CHECK( max( sqrAbs(psiAllv-psiAllvv) ) < 1e-24 );
  // For some reason, there is no exact equality


  ////////////////////////////////////////
  // Composite trafo with one non-Identity
  ////////////////////////////////////////

  trafoFunction(psi23,psi23v);

  {
    using namespace blitz::tensor;
    psiAllv=psi1(i)*psi23v(j,k)*psi4(l)*psi56(m,n);
  }

  {
    using namespace boost::fusion;

    transformation::Composite<result_of::make_vector<Identity1,void(*)(const Array2&, Array2&),Identity1,Identity2>::type>
      (make_vector(Identity1(),&trafoFunction,Identity1(),Identity2())).transform(psiAll,psiAllvv);
  }


  BOOST_CHECK( max( sqrAbs(psiAllv-psiAllvv) ) < 1e-24 );

}

