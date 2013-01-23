// -*- C++ -*-
#ifndef UTILS_INCLUDE_MULTIINDEXITERATOR_H_INCLUDED
#define UTILS_INCLUDE_MULTIINDEXITERATOR_H_INCLUDED

#include "MultiIndexIteratorFwd.h"

#include "BlitzTiny.h"

#include <boost/operators.hpp>

#include <boost/mpl/identity.hpp>


namespace cpputils {


namespace details {


template<int RANK>
struct IterHelper : boost::input_iterator_helper<MultiIndexIterator<RANK>,TTD_IDXTINY(RANK)> {};

template<int RANK>
void doIt(const TTD_IDXTINY(RANK)& lbound, 
	  const TTD_IDXTINY(RANK)& ubound,
	  TTD_IDXTINY(RANK)& idx,
	  boost::mpl::int_<0>);


template<int RANK, int N>
void doIt(const TTD_IDXTINY(RANK)& lbound, 
	  const TTD_IDXTINY(RANK)& ubound,
	  TTD_IDXTINY(RANK)& idx,
	  boost::mpl::int_<N>);


} // details


template<int RANK>
class MultiIndexIterator 
  : public details::IterHelper<RANK>
{
public:
  typedef details::IterHelper<RANK> Base;

  typedef typename Base::value_type IdxTiny;
  
  typedef boost::mpl::false_ Begin;
  typedef boost::mpl::true_  End  ;
  
  static const Begin begin;
  static const End   end  ;

  MultiIndexIterator(const IdxTiny& lbound, const IdxTiny& ubound, Begin)
    : lbound_(lbound), ubound_(ubound), idx_(lbound_) {}
  // if it's not the end, it's the beginning
  MultiIndexIterator(const IdxTiny& lbound, const IdxTiny& ubound, End  )
    : lbound_(lbound), ubound_(ubound), idx_(ubound_) {operator++();}
  // the end iterator is in fact beyond the end by one, hence the need for the increment
  
  MultiIndexIterator& operator=(const MultiIndexIterator& other) {idx_=other.idx_; return *this;}
  
  MultiIndexIterator& operator++() {details::doIt(lbound_,ubound_,idx_,boost::mpl::int_<RANK-1>()); return *this;}

  const IdxTiny& operator*() const {return idx_;}
        IdxTiny& operator*()       {return const_cast<IdxTiny&>(static_cast<const MultiIndexIterator*>(this)->operator*());}

  friend bool operator==(const MultiIndexIterator& i1, const MultiIndexIterator& i2) {return all(i1.idx_==i2.idx_);}
  // The user has to take care that the bounds are actually the same

  const MultiIndexIterator getBegin() const {return MultiIndexIterator(lbound_,ubound_,begin);}
  const MultiIndexIterator getEnd  () const {return MultiIndexIterator(lbound_,ubound_,end  );}
  
  MultiIndexIterator& setToBegin() {idx_=lbound_; return    *this ;}
  MultiIndexIterator& setToEnd  () {idx_=ubound_; return ++(*this);}
  
private:
  const IdxTiny lbound_, ubound_;

  IdxTiny idx_;

};


} // cpputils



#endif // UTILS_INCLUDE_MULTIINDEXITERATOR_H_INCLUDED
