// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFile{Defines class MultiIndexIterator and a few helpers}
#ifndef CPPQEDCORE_UTILS_MULTIINDEXITERATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_MULTIINDEXITERATOR_H_INCLUDED

#include "MultiIndexIteratorFwd.h"

#include "BlitzTiny.h"

#include <boost/operators.hpp>

#include <boost/mpl/identity.hpp>


/// Namespace comprising otherwise hard-to-classify generic utilities
namespace cpputils {


/// Helpers to MultiIndexIterator
namespace mii {
 
typedef boost::mpl::false_ Begin; 
typedef boost::mpl::true_  End  ;
  
const Begin begin=Begin();
const End   end  =End  ();
  
}


#define FORWARD_IteratorHelper boost::forward_iterator_helper<MultiIndexIterator<RANK>,IdxTiny<RANK> >

/// An iterator over all possible combinations of indices (a number of integers) between certain bounds.
/**
 * Models an [input iterator](http://www.cplusplus.com/reference/std/iterator/InputIterator/).
 * 
 * \tparam RANK the number of indices
 * 
 */
template<int RANK>
class MultiIndexIterator : public FORWARD_IteratorHelper
{
private:
  typedef FORWARD_IteratorHelper Base;
  
#undef FORWARD_IteratorHelper

public:
  typedef typename Base::value_type MultiIndex;
  typedef typename MultiIndex::T_numtype SingleIndex;

  /// \name Constructors
  //@{
    /// Initialization to the beginning of the sequence
  MultiIndexIterator(const MultiIndex& lbound, ///< tiny vector comprising the sequence of lower bounds
                     const MultiIndex& ubound, ///< tiny vector comprising the sequence of upper bounds (inclusive!)
                     mii::Begin)
    : lbound_(lbound), ubound_(ubound), idx_(lbound_) {}

  MultiIndexIterator(const MultiIndex& lbound, const MultiIndex& ubound, mii::End  )
    : lbound_(lbound), ubound_(ubound), idx_(ubound_) {operator++();} ///< ” to the end of the sequence (which is in fact beyond the end by one)

    /// \overload
    /** `lbound` is all-zero */
  MultiIndexIterator(const MultiIndex& ubound, ///< tiny vector comprising the sequence of upper bounds (inclusive!)
                     mii::Begin b)
    : MultiIndexIterator(MultiIndex(SingleIndex(0)),ubound,b) {}

    /// \overload
    /** `lbound` is all-zero */
  MultiIndexIterator(const MultiIndex& ubound, mii::End e)
    : MultiIndexIterator(MultiIndex(SingleIndex(0)),ubound,e) {}
  //@}
  
  MultiIndexIterator& operator=(const MultiIndexIterator& other) {idx_=other.idx_; return *this;}
  
  /// \name InputIterator operations
  //@{
  MultiIndexIterator& operator++() {doIt(boost::mpl::int_<RANK-1>()); return *this;}

  const MultiIndex& operator*() const {return idx_;}
        MultiIndex& operator*()       {return const_cast<MultiIndex&>(static_cast<const MultiIndexIterator*>(this)->operator*());}
  //@}

  /// Comparison only for the actual values referred to. The user has to take care that the bounds are actually the same!
  friend bool operator==(const MultiIndexIterator& i1, const MultiIndexIterator& i2) {return all(i1.idx_==i2.idx_);}

  /// \name Convenience
  //@{
  const MultiIndexIterator getBegin() const {return MultiIndexIterator(lbound_,ubound_,mii::begin);}
  const MultiIndexIterator getEnd  () const {return MultiIndexIterator(lbound_,ubound_,mii::end  );}
  
  MultiIndexIterator& setToBegin() {idx_=lbound_; return    *this ;}
  MultiIndexIterator& setToEnd  () {idx_=ubound_; return ++(*this);}
  //@}
  
private:
  void doIt(boost::mpl::int_<0>);

  template<int N>
  void doIt(boost::mpl::int_<N>);

  const MultiIndex lbound_, ubound_;

  MultiIndex idx_;

};


} // cpputils



#endif // CPPQEDCORE_UTILS_MULTIINDEXITERATOR_H_INCLUDED
