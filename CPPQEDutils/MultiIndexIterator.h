// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines class MultiIndexIterator and a few helpers}
#ifndef CPPQEDCORE_UTILS_MULTIINDEXITERATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_MULTIINDEXITERATOR_H_INCLUDED

#include "MultiIndexIterator.h"

#include "BlitzTiny.h"

#include <boost/operators.hpp>

#include <boost/mpl/identity.hpp>


/// Namespace comprising otherwise hard-to-classify generic utilities
namespace cppqedutils {


/// Helpers to MultiIndexIterator
namespace mii {
 
typedef std::false_type Begin; 
typedef std::true_type  End  ;
  
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

  /// \name Constructors
  //@{
  template <bool IS_END>
  MultiIndexIterator(const MultiIndex& lbound, ///< tiny vector comprising the sequence of lower bounds
                     const MultiIndex& ubound, ///< tiny vector comprising the sequence of upper bounds (inclusive!)
                     std::bool_constant<IS_END>)
    : lbound_(lbound), ubound_(ubound), idx_(IS_END ? ubound_ : lbound_) {if constexpr (IS_END) operator++();}

    /// \overload
    /** `lbound` is all-zero */
  template <bool IS_END>
  MultiIndexIterator(const MultiIndex& ubound, ///< tiny vector comprising the sequence of upper bounds (inclusive!)
                     std::bool_constant<IS_END> t)
    : MultiIndexIterator(MultiIndex(typename MultiIndex::T_numtype{0}),ubound,t) {}
  //@}
  
  MultiIndexIterator& operator=(const MultiIndexIterator&) = delete; // {idx_=other.idx_; return *this;}
  
  /// \name InputIterator operations
  //@{
  MultiIndexIterator& operator++() {doIt<RANK-1>(); return *this;}

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
  template<int N>
  void doIt()
  {
    if constexpr (N!=0) {
      if (idx_(N)==ubound_(N)) {idx_(N)=lbound_(N); doIt<N-1>();}
      else idx_(N)++;
    }
    else idx_(0)++; // This will of course eventually put the iterator into an illegal state when idx(0)>ubound(0), but this is how every (unchecked) iterator works.
  }

  const MultiIndex lbound_, ubound_;

  MultiIndex idx_;

};


} // cppqedutils



#endif // CPPQEDCORE_UTILS_MULTIINDEXITERATOR_H_INCLUDED
