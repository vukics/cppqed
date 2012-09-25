// -*- C++ -*-
#ifndef   UTILS_INCLUDE_SMARTPTR_H_INCLUDED
#define   UTILS_INCLUDE_SMARTPTR_H_INCLUDED

#include "SmartPtrFwd.h"

#include <boost/shared_ptr.hpp>


namespace cpputils {


class NullDeleter
{
public:
  void operator()(const void*) const {}

};


template<typename T>
const boost::shared_ptr<T>
nonOwningSharedPtr(T* t)
{
  return boost::shared_ptr<T>(t,NullDeleter());
}


template<typename T>
const boost::shared_ptr<const T>
nonOwningConstSharedPtr(T* t)
{
  return boost::shared_ptr<const T>(t,NullDeleter());
}


template<typename T>
boost::shared_ptr<T> sharedPointerize(boost::shared_ptr<T> t) {return                               t ;}

template<typename T>
boost::shared_ptr<T> sharedPointerize(                  T& t) {return cpputils::nonOwningSharedPtr(&t);}

template<typename T>
boost::shared_ptr<T> sharedPointerize(                  T* t) {return cpputils::nonOwningSharedPtr( t);}

/*
template<typename T>
boost::shared_ptr<const T> sharedPointerize(boost::shared_ptr<const T> t) {return                               t ;}

template<typename T>
boost::shared_ptr<const T> sharedPointerize(const                   T& t) {return cpputils::nonOwningSharedPtr(&t);}

template<typename T>
boost::shared_ptr<const T> sharedPointerize(const                   T* t) {return cpputils::nonOwningSharedPtr( t);}
*/

} // cpputils


#endif // UTILS_INCLUDE_SMARTPTR_H_INCLUDED
