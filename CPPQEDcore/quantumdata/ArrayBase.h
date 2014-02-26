// -*- C++ -*-
/// \briefFileDefault
#ifndef   CPPQEDCORE_QUANTUMDATA_ARRAYBASE_H_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_ARRAYBASE_H_INCLUDED

#include "ArrayBaseFwd.h"

#include "BlitzArrayExtensions.h"
#include "CMatrix.h"
#include "ComplexArrayExtensions.h"

#include <boost/utility.hpp>


namespace quantumdata {


/// Comprises the common functionalities of StateVector and DensityOperator.
template<int RANK>
class ArrayBase : private boost::noncopyable
{
protected:
  typedef CArray<RANK> ArrayLow; ///< The underlying storage
  typedef linalg::CVector CVector;

  explicit ArrayBase(const ArrayLow& arrayLow) : arrayLow_(arrayLow) {} ///< By-reference semantics (basically the copy of a `blitz::Array`). Apart from this, copying is not possible.

  virtual ~ArrayBase() {}

  /// Assignment with by-value semantics (like the assignment of a `blitz::Array`).
  /**
   * Default assignment synthetised by the compiler ditto.
   */ 
  ArrayBase& operator=(const ArrayLow& arrayLow ) {arrayLow_=arrayLow; return *this;}
  
  /// \name The underlying ArrayLow
  //@{
  const ArrayLow& getArray() const {return arrayLow_;}
        ArrayLow& getArray()       {return const_cast<ArrayLow&>(static_cast<const ArrayBase*>(this)->getArray());}
  //@}
  
  /// \name Naive vector-space operations
  //@{
  void operator+=(const ArrayBase& arrayBase) {arrayLow_+=arrayBase.arrayLow_;}
  void operator-=(const ArrayBase& arrayBase) {arrayLow_-=arrayBase.arrayLow_;}
  //@}
  
  /// \name Naive vector-space operations allowing also for mixed-mode arithmetic
  //@{
  template<typename OTHER>
  void operator*=(const OTHER& dc) {arrayLow_*=dc;}

  template<typename OTHER>
  void operator/=(const OTHER& dc) {arrayLow_/=dc;}
  //@}

  /// \name One-dimensional view of the underlying data
  //@{
    /// 1d view created on the fly via blitzplusplus::unaryArray.
    /**
      * \note This is meant to be used only if the underlying storage is contiguous, which may not be the case since the object may reference arrays of any layout.
      * In debug mode, a non-contiguous storage is detected by the implementing function blitzplusplus::unaryArray, and an exception of type 
      * blitzplusplus::NonContiguousStorageException is thrown.
      */
  const CVector vectorView() const {return blitzplusplus::unaryArray(arrayLow_);}
        CVector vectorView()       {return blitzplusplus::unaryArray(arrayLow_);} ///< ”
  // The usual technique of defining the non-const in terms of the const doesn't work here, because `CVector is not a pointer, reference, nor a pointer-to-data-member type'
  //@}

  /// The entrywise array norm
  /**
   * \f[\norm A _{\text{F}}=\sqrt{\sum_i\,\abs{A_i}^2},\f]
   * with \f$i\f$ running through all the multi-indices.
   */
  double frobeniusNorm() const {return sqrt(sum(blitzplusplus::sqrAbs(arrayLow_)));}

private:
  ArrayLow arrayLow_;

};


} // quantumdata



#endif // CPPQEDCORE_QUANTUMDATA_ARRAYBASE_H_INCLUDED
