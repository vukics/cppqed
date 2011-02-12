// -*- C++ -*-
#ifndef   DISPATCH_FREE_TYPE_INCLUDED
#define   DISPATCH_FREE_TYPE_INCLUDED

#include <boost/smart_ptr.hpp>

/*
template<typename FT>
const FT*const dispatchFreeType(const             FT* free) {return  free      ;}
*/

template<typename FT>
const FT*const dispatchFreeType(boost::shared_ptr<FT> free) {return  free.get();}


template<typename FT>
const FT*const dispatchFreeType(const             FT& free) {return &free      ;}


struct EmptyAveragingBaseForBinaryInteractions {};


#endif // DISPATCH_FREE_TYPE_INCLUDED
