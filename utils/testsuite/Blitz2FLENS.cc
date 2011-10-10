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


#include "Blitz2FLENS.h"

#include "MathExtensions.h"
#include "BlitzTiny.h"
#include "BlitzArray.h"
#include "ComplexArrayExtensions.h"
#include "Randomized.h"

#include "Range.h"

#include<flens/flens.h>

#include<boost/bind.hpp>


using namespace std;
using namespace randomized;
using namespace blitz2flens;
using namespace blitzplusplus;
using namespace mathutils;

using blitz::ColumnMajorArray;

const double epsilonCmp=1e-12;

const int RA=3;


typedef TTD_CARRAY(  RA) CAR ;
typedef TTD_DARRAY(  RA) DAR ;
typedef TTD_CARRAY(2*RA) CA2R;


typedef DenseVectorMF<dcomp >::type CDenseVector;
typedef DenseVectorMF<double>::type DDenseVector;
typedef GeMatrixMF<dcomp,RowMajor>::type GeMatrixRM;
typedef GeMatrixMF<dcomp,ColMajor>::type GeMatrixCM;

typedef HeMatrixMF<RowMajor>::type HeMatrixRM;

int main()
{
  Randomized::SmartPtr ran(MakerGSL()(1001));

  TTD_EXTTINY(RA) dims(6,4,5);

  // RowMajor, C/C++Array

  {
    CAR  vecBlitz(dims);
    CA2R matBlitz(concatenateTinies(dims,dims));
    CA2R vlBlitz(matBlitz.shape()), vrBlitz(matBlitz.shape());

    CDenseVector vecFLENS(blitz2flens::vector(vecBlitz));

    GeMatrixRM
      vlFLENS(matrix(vlBlitz,RowMajorTag())),
      vrFLENS(matrix(vrBlitz,RowMajorTag()));


    fillWithRandom(matBlitz,ran);
    
    {
      CA2R a(matBlitz.copy());
      GeMatrixRM aFLENS(matrix(a,RowMajorTag()));
      cerr<<"Entering ev routine ... "; assert(!ev(true,true,aFLENS,vecFLENS,vlFLENS,vrFLENS)); cerr<<"exiting, checking result ... ";
    }

    {
      CA2R resTensor(matBlitz.shape());
      TTD_CARRAY(9) temp9(concatenateTinies(matBlitz.shape(),dims));
      {
	using namespace blitz::tensor;
	temp9=matBlitz(i,j,k,o,p,q)*conj(vlBlitz(l,m,n,o,p,q))/vecBlitz(l,m,n);
	resTensor=CA2R(sum(sum(sum(temp9,q),p),o));

	assert(!fcmp(1-max(abs(vlBlitz.transpose(3,4,5,0,1,2)-conj(resTensor))),1,epsilonCmp));

	temp9=matBlitz(o,p,q,i,j,k)*vrBlitz(l,m,n,o,p,q)/vecBlitz(l,m,n);
	resTensor=CA2R(sum(sum(sum(temp9,q),p),o));
      }
      assert(!fcmp(1-max(abs(vrBlitz.transpose(3,4,5,0,1,2)-resTensor)),1,epsilonCmp));

      cerr<<"Nonsymmetric eigenproblem in RowMajor OK!\n";

    }
  }

  // ColMajor, FortranArray

  {
    CAR  vecBlitz(dims,ColumnMajorArray<RA>());
    CA2R matBlitz(concatenateTinies(dims,dims),ColumnMajorArray<2*RA>());
    CA2R vlBlitz(matBlitz.shape(),ColumnMajorArray<2*RA>()), vrBlitz(matBlitz.shape(),ColumnMajorArray<2*RA>());

    CDenseVector vecFLENS(blitz2flens::vector(vecBlitz));

    GeMatrixCM
      vlFLENS(matrix(vlBlitz,ColMajorTag())),
      vrFLENS(matrix(vrBlitz,ColMajorTag()));


    fillWithRandom(matBlitz,ran);
    
    {
      CA2R a(matBlitz.copy());
      GeMatrixCM aFLENS(matrix(a,ColMajorTag()));
    
      cerr<<"Entering ev routine ... "; assert(!ev(true,true,aFLENS,vecFLENS,vlFLENS,vrFLENS)); cerr<<"exiting, checking result ... ";
    }

    {
      CA2R resTensor(matBlitz.shape());
      TTD_CARRAY(9) temp9(concatenateTinies(matBlitz.shape(),dims));
      {
	using namespace blitz::tensor;
	temp9=matBlitz(i,j,k,o,p,q)*vrBlitz(o,p,q,l,m,n)/vecBlitz(l,m,n);
	resTensor=CA2R(sum(sum(sum(temp9,q),p),o));

	assert(!fcmp(1-max(abs(vrBlitz-resTensor)),1,epsilonCmp));

	temp9=matBlitz(o,p,q,i,j,k)*conj(vlBlitz(o,p,q,l,m,n))/vecBlitz(l,m,n);
	resTensor=CA2R(sum(sum(sum(temp9,q),p),o));
      }

      assert(!fcmp(1-max(abs(vlBlitz-conj(resTensor))),1,epsilonCmp));

      cerr<<"\"                            ColMajor OK!\n";

    }
  }

  {
    CAR  vecBlitz1(dims);
    DAR  vecBlitz2(dims);
    CA2R matBlitz(concatenateTinies(dims,dims));

    fillWithRandom(matBlitz,ran);

    matBlitz+=hermitianConjugate(matBlitz);
    
    {
      CA2R a(matBlitz.copy());

      CDenseVector vecFLENS(blitz2flens::vector(vecBlitz1));
      GeMatrixRM aFLENS(matrix(a,RowMajorTag()));
    
      cerr<<"Entering ev routine ... "; ev(false,false,aFLENS,vecFLENS,aFLENS,aFLENS); cerr<<"exiting, checking result ... ";
    }
    assert(!fcmp(1-max(abs(imag(vecBlitz1))),1,epsilonCmp));
    cerr<<"Hermitian eigenproblem eigenvalues all real.\n";

    {
      CA2R a(matBlitz.copy());

      DDenseVector vecFLENS(blitz2flens::vector(vecBlitz2));
      HeMatrixRM aFLENS(hermitianMatrix(a,RowMajorTag()));
    
      cerr<<"Entering ev routine ... "; ev(true,aFLENS,vecFLENS); cerr<<"exiting, checking result ... ";

      {
	CA2R resTensor(matBlitz.shape());
	TTD_CARRAY(9) temp9(concatenateTinies(matBlitz.shape(),dims));
	{
	  using namespace blitz::tensor;
	  temp9=matBlitz(o,p,q,i,j,k)*a(l,m,n,o,p,q)/vecBlitz2(l,m,n);
	  resTensor=CA2R(sum(sum(sum(temp9,q),p),o));
	  
	  assert(!fcmp(1-max(abs(a.transpose(3,4,5,0,1,2)-resTensor)),1,epsilonCmp));
	  
	}

      }

    }

    cerr<<"\"                      eigenproblem OK.\n";

    std::sort(vecBlitz1.data(),vecBlitz1.data()+vecBlitz1.size(),realCompare);

    // note: vecBlitz2 is already sorted

    assert(!fcmp(1-max(abs(vecBlitz2-real(vecBlitz1))),1,epsilonCmp));

    cerr<<"Hermitian eigenproblem eigenvalues match.\n";

  }

}

