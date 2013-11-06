// -*- C++ -*-
#ifndef UTILS_FUNCTIONAL_H_INCLUDED
#define UTILS_FUNCTIONAL_H_INCLUDED

#include<functional>


namespace cpputils {


template<typename T>
struct plus : public std::binary_function<T,T,T>
{
  T operator()(T t1, T t2) const {return T(t1+t2);}
};


} // cpputils


#endif // UTILS_FUNCTIONAL_H_INCLUDED
