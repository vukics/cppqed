// -*- C++ -*-
#ifndef _CPPUTILS_FUNCTIONAL_H
#define _CPPUTILS_FUNCTIONAL_H

#include<functional>


namespace cpputils {


template<typename T>
struct plus : public std::binary_function<T,T,T>
{
  T operator()(T t1, T t2) const {return T(t1+t2);}
};


} // cpputils


#endif // _CPPUTILS_FUNCTIONAL_H
