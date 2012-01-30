// -*- C++ -*-
#ifndef   ARRAY_BASE_INCLUDED
#define   ARRAY_BASE_INCLUDED

#include "ArrayBaseFwd.h"

#include "CMatrix.h"
#include "BlitzArrayExtensions.h"


namespace quantumdata {



template<int RANK>
class ArrayBase : private boost::noncopyable
{
protected:
  typedef TTD_CARRAY(RANK) ArrayLow;
  typedef linalg::CVector CVector;

  explicit ArrayBase(const ArrayLow& arrayLow) : arrayLow_(arrayLow) {} // By reference semantics

  virtual ~ArrayBase() {}

  // Assignment: by value semantics. 
  ArrayBase& operator=(const ArrayLow& arrayLow ) {arrayLow_=arrayLow; return *this;}
  // + default assignment

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

  const CVector vectorView() const {return blitzplusplus::unaryArray(arrayLow_);}
  CVector       vectorView()       {return blitzplusplus::unaryArray(arrayLow_);}
  // The usual technique doesn't work here, because `CVector is not a pointer, reference, nor a pointer-to-data-member type'


private:
  ArrayLow arrayLow_;

};


} // quantumdata



#endif // ARRAY_BASE_INCLUDED
