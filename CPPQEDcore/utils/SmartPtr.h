// -*- C++ -*-
#ifndef   UTILS_SMARTPTR_H_INCLUDED
#define   UTILS_SMARTPTR_H_INCLUDED

#include <boost/shared_ptr.hpp>


namespace cpputils {

template<typename T>
inline void nullDelete (T*) {}


template<typename T>
const boost::shared_ptr<T>
nonOwningSharedPtr(T* t)
{
  return boost::shared_ptr<T>(t,nullDelete<T>);
}


template<typename T>
const boost::shared_ptr<const T>
nonOwningConstSharedPtr(T* t)
{
  return boost::shared_ptr<const T>(t,nullDelete<T>);
}


template<typename T>
const boost::shared_ptr<T> sharedPointerize(boost::shared_ptr<T> t) {return                     t ;}

template<typename T>
const boost::shared_ptr<T> sharedPointerize(                  T& t) {return nonOwningSharedPtr(&t);}

template<typename T>
const boost::shared_ptr<T> sharedPointerize(                  T* t) {return nonOwningSharedPtr( t);}


} // cpputils


#endif // UTILS_SMARTPTR_H_INCLUDED
