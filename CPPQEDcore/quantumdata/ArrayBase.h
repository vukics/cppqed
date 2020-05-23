// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_QUANTUMDATA_ARRAYBASE_H_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_ARRAYBASE_H_INCLUDED

#include "BlitzArrayExtensions.h"
#include "CMatrix.h"
#include "ComplexArrayExtensions.h"
#include "Operators.h"

#include <boost/utility.hpp>


namespace quantumdata {

struct ByReference {}; const ByReference byReference{};

template<typename>
struct ArrayRank;

/// Comprises the common functionalities of StateVector and DensityOperator.
template<typename Derived>
class ArrayBase : private boost::noncopyable, private linalg::VectorSpace<Derived>
{
protected:
  typedef CArray<ArrayRank<Derived>::value> ArrayLow; ///< The underlying storage
  typedef linalg::CVector CVector;

  ArrayBase() : arrayLow_(0) {}
  
  explicit ArrayBase(const ArrayLow& arrayLow) : arrayLow_(arrayLow) {} ///< By-reference semantics (basically the copy of a `blitz::Array`). Apart from this, copying is not possible.

  ArrayBase(ArrayLow&& array) : arrayLow_(array) {}
  
  ArrayBase(ArrayBase&& array) : ArrayBase(std::move(array.getArray())) {}
  
  virtual ~ArrayBase() {}

public:
  /// Mixed-mode assignment with by-value semantics
  /**
  * The standard assignment and the templated assignment together cover a lot of possibilities, including also assignment from a StateVectorLow,
  * but for example also from a DArray<RANK>, or just a const c-number. (Can be assigned from anything a CArray<RANK> can be assigned from.)
  * 
  * \tparam OTHER the “other” type in mixed mode
  */
  template <typename OTHER>
  Derived& operator=(const OTHER& other) {arrayLow_=other; return static_cast<Derived&>(*this);}
  
  /// \name The underlying ArrayLow
  //@{
  const ArrayLow& getArray() const {return arrayLow_;}
        ArrayLow& getArray()       {return const_cast<ArrayLow&>(static_cast<const ArrayBase*>(this)->getArray());}
  //@}
  
  /// \name Naive vector-space operations
  //@{
  Derived& operator+=(const ArrayBase& arrayBase) {arrayLow_+=arrayBase.arrayLow_; return static_cast<Derived&>(*this);}
  Derived& operator-=(const ArrayBase& arrayBase) {arrayLow_-=arrayBase.arrayLow_; return static_cast<Derived&>(*this);}

  Derived operator-() const {Derived res(this->getDimensions(),false); res.getArray()=-this->getArray(); return res;} ///< involves a deep-copy
  Derived operator+() const {return *this;} ///< simply deep copy
  //@}

  /// \name Naive vector-space operations allowing also for mixed-mode arithmetics
  //@{
  template<typename OTHER>
  Derived& operator*=(const OTHER& dc) {arrayLow_*=dc; return static_cast<Derived&>(*this);} ///< \tparam OTHER the “other” type in mixed mode

  template<typename OTHER>
  Derived& operator/=(const OTHER& dc) {arrayLow_/=dc; return static_cast<Derived&>(*this);}
  //@}

protected:
  /// \name One-dimensional view of the underlying data
  //@{
    /// 1d view created on the fly via blitzplusplus::unaryArray.
    /**
      * \note This is meant to be used only if the underlying storage is contiguous, which may not be the case since the object may reference arrays of any layout.
      * In debug mode, a non-contiguous storage is detected by the implementing function blitzplusplus::unaryArray, and an exception of type 
      * blitzplusplus::NonContiguousStorageException is thrown.
      */
  auto vectorView() const {return blitzplusplus::unaryArray(arrayLow_);}
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
