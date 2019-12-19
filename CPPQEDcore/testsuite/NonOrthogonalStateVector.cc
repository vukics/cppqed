// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#define BOOST_TEST_MODULE StateVector
#include <boost/test/unit_test.hpp>

#include "NonOrthogonalStateVector.h"

using namespace quantumdata;

typedef CArray<1> Array1;
typedef CArray<2> Array2;
typedef NonOrthogonalStateVector<1, Array2> NOSV_Array2;
typedef transformation::Identity<1> I1;
typedef NonOrthogonalStateVector<1, I1> NOSV_I1;
typedef transformation::Identity<2> I2;
typedef NonOrthogonalStateVector<2, I2> NOSV_I2;

BOOST_AUTO_TEST_CASE(Dual)
{
    // Create Transformation Matrix
    Array2 T(3);
    T = 1,0,0,
        0,1,dcomp(0.5,-0.3),
        0,dcomp(0.5,0.3),1;

    // Create StateVector
    Array1 x(3);
    x = 1,2,3;
    NOSV_Array2 sv1(x,T,byReference);

    // Calculate and check dual vector
    sv1.update();
    Array1 dual = sv1.dual();
    BOOST_CHECK_EQUAL(dual(0), x(0)*T(0,0)+x(1)*T(0,1)+x(2)*T(0,2));
    BOOST_CHECK_EQUAL(dual(1), x(0)*T(1,0)+x(1)*T(1,1)+x(2)*T(1,2));
    BOOST_CHECK_EQUAL(dual(2), x(0)*T(2,0)+x(1)*T(2,1)+x(2)*T(2,2));

    // Check norm
    BOOST_CHECK_EQUAL(sv1.norm(), sqrt(real(
                                    sv1()(0)*conj(dual(0)) +
                                    sv1()(1)*conj(dual(1)) +
                                    sv1()(2)*conj(dual(2))
                                    )));
};

BOOST_AUTO_TEST_CASE(CompositeStateVector)
{
    // Create Transformation Matrix
    Array2 T(3,3);
    T = 1,0,0,
        0,1,dcomp(0.5,-0.5),
        0,dcomp(0.5,0.5),1;

    // Create StateVectors
    Array1 x(3);
    x = 1,2,3;
    NOSV_Array2 sv1(x,T,byReference);

    Array2 y(3,3);
    y = 1,2,
        3,1;
    I2 i2;
    NOSV_I2 sv2(y,i2,byReference);

    typedef TensorType<NOSV_I2, NOSV_Array2>::type Result;
    typedef Result::StateVectorLow ResultArray;
    Result res1(sv2*sv1);

    // Calculate all dual vectors.
    res1.update();
    sv1.update();
    sv2.update();

    // Check Composite StateVector
    ResultArray res2(blitzplusplus::doDirect(sv2.dual(), sv1.dual(),
                                    blitzplusplus::dodirect::Mul()));
    blitz::Array<bool,3> b(3,3,3);
    BOOST_CHECK(blitz::all(res1.dual()==res2));
}

BOOST_AUTO_TEST_CASE(vectorspace_operations)
{
    Array1 x(3), y(3);
    x=1,2,3;
    y=2,3,1;
    I1 i1;
    NOSV_I1 sv1(x, i1, byReference);
    NOSV_I1 sv2(y, i1, byReference);
    NOSV_I1 result_p(sv1 + sv2);
    BOOST_CHECK_EQUAL(result_p()(0), x(0)+y(0));
    BOOST_CHECK_EQUAL(result_p()(1), x(1)+y(1));
    BOOST_CHECK_EQUAL(result_p()(2), x(2)+y(2));

    NOSV_I1 result_m(sv1 - sv2);
    BOOST_CHECK_EQUAL(result_m()(0), x(0)-y(0));
    BOOST_CHECK_EQUAL(result_m()(1), x(1)-y(1));
    BOOST_CHECK_EQUAL(result_m()(2), x(2)-y(2));

    NOSV_I1 result_mul(sv1*dcomp(4));
    BOOST_CHECK_EQUAL(result_mul()(0), x(0)*dcomp(4));
    BOOST_CHECK_EQUAL(result_mul()(1), x(1)*dcomp(4));
    BOOST_CHECK_EQUAL(result_mul()(2), x(2)*dcomp(4));
}



BOOST_AUTO_TEST_CASE(tensorproduct)
{
    I1 i1;
    I2 i2;
    Array1 x(3);
    x = 1,2,3;
    Array2 y(3,3);
    y = 4,3,2,
        5,4,3,
        6,1,2;

    NOSV_I1 sv1(x, i1, byReference);
    NOSV_I2 sv2(y, i2, byReference);

    typedef TensorType<NOSV_I1,NOSV_I2>::type NOSV_I1I2;
    NOSV_I1I2 sv_i1i2(sv1*sv2);
    typedef TensorType<NOSV_I2,NOSV_I1>::type NOSV_I2I1;
    NOSV_I2I1 sv_i2i1(sv2*sv1);
    typedef TensorType<NOSV_I1I2,NOSV_I1>::type NOSV_I1I2I1_a;
    typedef TensorType<NOSV_I1,NOSV_I2I1>::type NOSV_I1I2I1_b;
    BOOST_MPL_ASSERT(( boost::is_same<NOSV_I1I2I1_a, NOSV_I1I2I1_b> ));
    typedef NOSV_I1I2I1_a NOSV_I1I2I1;
    NOSV_I1I2I1(sv_i1i2*sv1);
    NOSV_I1I2I1(sv1*sv_i2i1);

    typedef StateVector<1> OSV;
    OSV osv(x, ByReference());
    NOSV_I2I1 mixed_a(sv2*osv);
    NOSV_I1I2 mixed_b(osv*sv2);
    typedef TensorType<NOSV_I2, OSV>::type NOSV_I2I1_orth;
    typedef TensorType<OSV, NOSV_I2>::type NOSV_I1I2_orth;
}



