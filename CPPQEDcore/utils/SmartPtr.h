// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Tools for creating non-owning shared pointers}
#ifndef   CPPQEDCORE_UTILS_SMARTPTR_H_INCLUDED
#define   CPPQEDCORE_UTILS_SMARTPTR_H_INCLUDED

#include <boost/shared_ptr.hpp>


namespace cpputils {


namespace details {

template<typename T>
inline void nullDelete (T*) {}

} // details


/// Returns a “non-owning” shared pointer that doesn’t do anything on deletion
/**
 * Can be used even on automatic variables. 
 * 
 * \tparam T the type of the variable
 * 
 * \warning Handle with due care! The wrapped object must live in a larger scope and must be properly deleted after the shared wrapper runs out of scope.
 * Failure in the first respect will result in dangling pointers, while the second in memory leak.
 * 
 */
template<typename T>
const boost::shared_ptr<T>
nonOwningSharedPtr(T* t ///< the object to be wrapped into the shared-pointer interface
                  )
{
  return boost::shared_ptr<T>(t,details::nullDelete<T>);
}


/// Returns a “non-owning” shared pointer to const that doesn’t do anything on deletion \copydetails nonOwningSharedPtr
template<typename T>
const boost::shared_ptr<const T>
nonOwningConstSharedPtr(T* t ///< the object to be wrapped into the shared-pointer-to-const interface
                       )
{
  return boost::shared_ptr<const T>(t,details::nullDelete<T>);
}


/// Part of a bundle of functions providing a unified interface to wrap objects into the shared-pointer interface, it simply returns its argument. \todo How to shared-pointerize an rvalue reference?
template<typename T>
const boost::shared_ptr<T> sharedPointerize(boost::shared_ptr<T> t) {return                     t ;}

/// Part of a bundle of functions providing a unified interface to wrap objects into the shared-pointer interface, it returns a \link nonOwningSharedPtr non-owning shared pointer\endlink to its argument.
template<typename T>
const boost::shared_ptr<T> sharedPointerize(                  T& t) {return nonOwningSharedPtr(&t);}

/// Part of a bundle of functions providing a unified interface to wrap objects into the shared-pointer interface, it returns a \link nonOwningSharedPtr non-owning shared pointer\endlink to its argument.
template<typename T>
const boost::shared_ptr<T> sharedPointerize(                  T* t) {return nonOwningSharedPtr( t);}


} // cpputils


#endif // CPPQEDCORE_UTILS_SMARTPTR_H_INCLUDED
