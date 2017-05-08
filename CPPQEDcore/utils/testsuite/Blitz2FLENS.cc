// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/*

The example demonstrates the use of LAPACK (via FLENS) with different
storage orders. Note that column major ordering of A is equivalent to
row major ordering of A^T.

right eigen:

A v = lambda v

left eigen:

u^H A = lambda u^H


      | A | A^T
----------------
left  | u | v^* 
----------------
right | v | u^* 
----------------

*/


#include "Blitz2FLENS.tcc"

#include "MathExtensions.h"
#include "BlitzTiny.h"
#include "ComplexArrayExtensions.tcc"
#include "Randomized.h"

#include "Range.h"

#include <flens/flens.h>

#define BOOST_TEST_MODULE Blitz2FLENS test
#include <boost/test/unit_test.hpp>

#include <boost/bind.hpp>


using namespace std;
using namespace randomized;
using namespace blitz2flens;
using namespace blitzplusplus;
using namespace mathutils;

using blitz::ColumnMajorArray;

const double epsilonCmp=1e-12;

const int RANK=3;


typedef CArray<  RANK> CAR ;
typedef DArray<  RANK> DAR ;
typedef CArray<2*RANK> CA2R;


typedef DenseVectorMF<dcomp >::type CDenseVector;
typedef DenseVectorMF<double>::type DDenseVector;
typedef GeMatrixMF<dcomp,RowMajor>::type GeMatrixRM;
typedef GeMatrixMF<dcomp,ColMajor>::type GeMatrixCM;

typedef HeMatrixMF<RowMajor>::type HeMatrixRM;


Randomized::Ptr ran(MakerGSL()(1001));

ExtTiny<RANK> dims(6,4,5);


// RowMajor, C/C++Array

BOOST_AUTO_TEST_CASE( RowMajorTest )
{
  CAR  vecBlitz(dims);
  CA2R matBlitz(concatenateTinies(dims,dims));
  CA2R vlBlitz(matBlitz.shape()), vrBlitz(matBlitz.shape());

  CDenseVector vecFLENS(blitz2flens::vector(vecBlitz));

  GeMatrixRM
    vlFLENS(matrix<RowMajor>(vlBlitz)),
    vrFLENS(matrix<RowMajor>(vrBlitz));


  fillWithRandom(matBlitz,ran);

  {
    CA2R a(matBlitz.copy());
    GeMatrixRM aFLENS(matrix<RowMajor>(a));
    cerr<<"Entering ev routine ... "; BOOST_CHECK(!ev(true,true,aFLENS,vecFLENS,vlFLENS,vrFLENS)); cerr<<"exiting, checking result ... ";
  }

  {
    CA2R resTensor(matBlitz.shape());
    CArray<9> temp9(concatenateTinies(matBlitz.shape(),dims));
    {
      using namespace blitz::tensor;
      temp9=matBlitz(i,j,k,o,p,q)*conj(vlBlitz(l,m,n,o,p,q))/vecBlitz(l,m,n);
      resTensor=CA2R(sum(sum(sum(temp9,q),p),o));

      BOOST_CHECK(!fcmp(1-max(abs(vlBlitz.transpose(3,4,5,0,1,2)-conj(resTensor))),1,epsilonCmp));

      temp9=matBlitz(o,p,q,i,j,k)*vrBlitz(l,m,n,o,p,q)/vecBlitz(l,m,n);
      resTensor=CA2R(sum(sum(sum(temp9,q),p),o));
    }
    BOOST_CHECK(!fcmp(1-max(abs(vrBlitz.transpose(3,4,5,0,1,2)-resTensor)),1,epsilonCmp));

    cerr<<"Nonsymmetric eigenproblem in RowMajor OK!\n";

  }
}


// ColMajor, FortranArray

BOOST_AUTO_TEST_CASE( ColMajorTest )
{
  CAR  vecBlitz(dims,ColumnMajorArray<RANK>());
  CA2R matBlitz(concatenateTinies(dims,dims),ColumnMajorArray<2*RANK>());
  CA2R vlBlitz(matBlitz.shape(),ColumnMajorArray<2*RANK>()), vrBlitz(matBlitz.shape(),ColumnMajorArray<2*RANK>());

  CDenseVector vecFLENS(blitz2flens::vector(vecBlitz));

  GeMatrixCM
    vlFLENS(matrix<ColMajor>(vlBlitz)),
    vrFLENS(matrix<ColMajor>(vrBlitz));


  fillWithRandom(matBlitz,ran);
    
  {
    CA2R a(matBlitz.copy());
    GeMatrixCM aFLENS(matrix<ColMajor>(a));
    
    cerr<<"Entering ev routine ... "; BOOST_CHECK(!ev(true,true,aFLENS,vecFLENS,vlFLENS,vrFLENS)); cerr<<"exiting, checking result ... ";
  }

  {
    CA2R resTensor(matBlitz.shape());
    CArray<9> temp9(concatenateTinies(matBlitz.shape(),dims));
    {
      using namespace blitz::tensor;
      temp9=matBlitz(i,j,k,o,p,q)*vrBlitz(o,p,q,l,m,n)/vecBlitz(l,m,n);
      resTensor=CA2R(sum(sum(sum(temp9,q),p),o));

      BOOST_CHECK(!fcmp(1-max(abs(vrBlitz-resTensor)),1,epsilonCmp));

      temp9=matBlitz(o,p,q,i,j,k)*conj(vlBlitz(o,p,q,l,m,n))/vecBlitz(l,m,n);
      resTensor=CA2R(sum(sum(sum(temp9,q),p),o));
    }

    BOOST_CHECK(!fcmp(1-max(abs(vlBlitz-conj(resTensor))),1,epsilonCmp));

    cerr<<"\"                            ColMajor OK!\n";

  }
}



BOOST_AUTO_TEST_CASE( HermitianTest )
{
  CAR  vecBlitz1(dims);
  DAR  vecBlitz2(dims);
  CA2R matBlitz(concatenateTinies(dims,dims));

  fillWithRandom(matBlitz,ran);

  matBlitz+=hermitianConjugate(matBlitz);
    
  {
    CA2R a(matBlitz.copy());

    CDenseVector vecFLENS(blitz2flens::vector(vecBlitz1));
    GeMatrixRM aFLENS(matrix<RowMajor>(a));
    
    cerr<<"Entering ev routine ... "; ev(false,false,aFLENS,vecFLENS,aFLENS,aFLENS); cerr<<"exiting, checking result ... ";
  }
  BOOST_CHECK(!fcmp(1-max(abs(imag(vecBlitz1))),1,epsilonCmp));
  cerr<<"Hermitian eigenproblem eigenvalues all real.\n";

  {
    CA2R a(matBlitz.copy());

    DDenseVector vecFLENS(blitz2flens::vector(vecBlitz2));
    HeMatrixRM aFLENS(hermitianMatrix<RowMajor>(a));
    
    cerr<<"Entering ev routine ... "; ev(true,aFLENS,vecFLENS); cerr<<"exiting, checking result ... ";

    {
      CA2R resTensor(matBlitz.shape());
      CArray<9> temp9(concatenateTinies(matBlitz.shape(),dims));
      {
        using namespace blitz::tensor;
        temp9=matBlitz(o,p,q,i,j,k)*a(l,m,n,o,p,q)/vecBlitz2(l,m,n);
        resTensor=CA2R(sum(sum(sum(temp9,q),p),o));
          
        BOOST_CHECK(!fcmp(1-max(abs(a.transpose(3,4,5,0,1,2)-resTensor)),1,epsilonCmp));
          
      }

    }

  }

  cerr<<"\"                      eigenproblem OK.\n";

  std::sort(vecBlitz1.data(),vecBlitz1.data()+vecBlitz1.size(),realCompare);

  // note: vecBlitz2 is already sorted

  BOOST_CHECK(!fcmp(1-max(abs(vecBlitz2-real(vecBlitz1))),1,epsilonCmp));

  cerr<<"Hermitian eigenproblem eigenvalues match.\n";

}

