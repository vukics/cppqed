// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_DETAILS_DISPATCHFREETYPE_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_DETAILS_DISPATCHFREETYPE_H_INCLUDED

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


#endif // ELEMENTS_INTERACTIONS_DETAILS_DISPATCHFREETYPE_H_INCLUDED
