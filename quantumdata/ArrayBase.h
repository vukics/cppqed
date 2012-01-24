// -*- C++ -*-
#ifndef   ARRAY_BASE_INCLUDED
#define   ARRAY_BASE_INCLUDED

#include "ArrayBaseFwd.h"

#include "BlitzArrayExtensions.h"
#include "CMatrix.h"
#include "ComplexArrayExtensions.h"


namespace quantumdata {



template<int RANK>
class ArrayBase
{
protected:
  typedef TTD_CARRAY(RANK) ArrayLow;
  typedef linalg::CVector CVector;

  explicit ArrayBase(const ArrayLow& arrayLow) : arrayLow_(arrayLow) {} // By reference semantics

  virtual ~ArrayBase() {}

  // Assignment: by value semantics. 
  // (Default assignment is wasteful as it assigns the vectorView_ as well, even though it refers to the same data.)
  ArrayBase& operator=(const ArrayLow & arrayLow ) {arrayLow_=arrayLow; return *this;}
  ArrayBase& operator=(const ArrayBase& arrayBase) {return operator=(arrayBase());}

  const ArrayLow& operator()() const {return arrayLow_;}
  ArrayLow& operator()() {return const_cast<ArrayLow&>(static_cast<const ArrayBase&>(*this)());}

  // naive operations for vector space

  void operator+=(const ArrayBase& arrayBase) {arrayLow_+=arrayBase();}
  void operator-=(const ArrayBase& arrayBase) {arrayLow_-=arrayBase();}


  template<typename OTHER>
  void operator*=(const OTHER& dc) {arrayLow_*=dc;}

  template<typename OTHER>
  void operator/=(const OTHER& dc) {arrayLow_/=dc;}


  // The following two can be called only if the underlying storage is contigous

  const CVector vectorView() const {return blitzplusplus::rankOneArray(arrayLow_);}
  CVector       vectorView()       {return blitzplusplus::rankOneArray(arrayLow_);}
  // The usual technique doesn't work here, because `CVector is not a pointer, reference, nor a pointer-to-data-member type'

  double frobeniusNorm() const {return sqrt(sum(blitzplusplus::sqrAbs(arrayLow_)));}

private:
  ArrayBase(const ArrayBase&); // forbid copying
  
  ArrayLow arrayLow_;

};


} // quantumdata



#endif // ARRAY_BASE_INCLUDED
