// -*- C++ -*-
#ifndef UTILS_INCLUDE_MULTIINDEXITERATOR_H_INCLUDED
#define UTILS_INCLUDE_MULTIINDEXITERATOR_H_INCLUDED

#include "MultiIndexIteratorFwd.h"

#include "BlitzTiny.h"

#include <boost/operators.hpp>

#include <boost/mpl/identity.hpp>


namespace cpputils {


namespace mii {
 
typedef boost::mpl::false_ Begin;
typedef boost::mpl::true_  End  ;
  
const Begin begin=Begin();
const End   end  =End  ();
  
}


#define INPUT_IteratorHelper boost::input_iterator_helper<MultiIndexIterator<RANK>,TTD_IDXTINY(RANK)>


template<int RANK>
class MultiIndexIterator : public INPUT_IteratorHelper
{
public:
  typedef INPUT_IteratorHelper Base;
  
#undef INPUT_IteratorHelper

  typedef typename Base::value_type IdxTiny;
  
  MultiIndexIterator(const IdxTiny& lbound, const IdxTiny& ubound, mii::Begin)
    : lbound_(lbound), ubound_(ubound), idx_(lbound_) {}
  // if it's not the end, it's the beginning
  MultiIndexIterator(const IdxTiny& lbound, const IdxTiny& ubound, mii::End  )
    : lbound_(lbound), ubound_(ubound), idx_(ubound_) {operator++();}
  // the end iterator is in fact beyond the end by one, hence the need for the increment
  
  MultiIndexIterator& operator=(const MultiIndexIterator& other) {idx_=other.idx_; return *this;}
  
  MultiIndexIterator& operator++() {doIt(boost::mpl::int_<RANK-1>()); return *this;}

  const IdxTiny& operator*() const {return idx_;}
        IdxTiny& operator*()       {return const_cast<IdxTiny&>(static_cast<const MultiIndexIterator*>(this)->operator*());}

  friend bool operator==(const MultiIndexIterator& i1, const MultiIndexIterator& i2) {return all(i1.idx_==i2.idx_);}
  // The user has to take care that the bounds are actually the same

  const MultiIndexIterator getBegin() const {return MultiIndexIterator(lbound_,ubound_,mii::begin);}
  const MultiIndexIterator getEnd  () const {return MultiIndexIterator(lbound_,ubound_,mii::end  );}
  
  MultiIndexIterator& setToBegin() {idx_=lbound_; return    *this ;}
  MultiIndexIterator& setToEnd  () {idx_=ubound_; return ++(*this);}
  
private:
  void doIt(boost::mpl::int_<0>);

  template<int N>
  void doIt(boost::mpl::int_<N>);

  const IdxTiny lbound_, ubound_;

  IdxTiny idx_;

};


} // cpputils



#endif // UTILS_INCLUDE_MULTIINDEXITERATOR_H_INCLUDED