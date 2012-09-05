// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_DETAILS_DISPATCHFREETYPE_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_DETAILS_DISPATCHFREETYPE_H_INCLUDED

#include "SmartPtr.h"

/*
template<typename FT>
const FT*const dispatchFreeType(const             FT* free) {return  free      ;}
*/

template<typename FT>
boost::shared_ptr<const FT> dispatchFreeType(boost::shared_ptr<const FT> free) {return                               free ;}


template<typename FT>
boost::shared_ptr<const FT> dispatchFreeType(const                   FT& free) {return cpputils::nonOwningSharedPtr(&free);}


struct EmptyAveragingBaseForBinaryInteractions {};


#endif // ELEMENTS_INTERACTIONS_DETAILS_DISPATCHFREETYPE_H_INCLUDED
