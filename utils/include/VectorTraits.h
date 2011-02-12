// -*- C++ -*-
/*
http://www.josuttis.com/libbook/memory/myalloc.hpp.html
*/
#ifndef _VECTOR_TRAITS_H
#define _VECTOR_TRAITS_H

#include<vector>
#include<limits>
#include<iostream>
#include<complex>



template<class> struct VectorTraits;

template<>
struct VectorTraits<std::vector<double> >
{
  typedef std::vector<double,blablabla::allocator<double> > dVec;
  
  static size_t size(const dVec& v) {return v.size();}

  static const double* data(const dVec& v) {return v.size() ? &v[0] : 0;}
  static       double* data(      dVec& v) {return const_cast<double*>(data(static_cast<const dVec&>(v)));}

  static dVec create(double* y, const dVec& v) 
  {
    const double dummy=0;
    blablabla::allocator<double> alloc(y);
    std::vector<double,blablabla::allocator<double> > x(v.size(),dummy,alloc);     
    return x;
  }
  static const dVec create(const double* y, const dVec& v) {return create(const_cast<double*>(y),v);}
};

template<>
struct VectorTraits<std::vector<dcomp> >
{
  typedef std::vector<dcomp,blablabla::allocator<dcomp> > cVec;
  
  static size_t size(const cVec& v) {return (v.size()<<1);} //The size of the underlying double* storage

  static const double* data(const cVec& v) {return v.size() ? &real(v[0]) : 0;}
  static       double* data(      cVec& v) {return const_cast<double*>(data(static_cast<const cVec&>(v)));}

  static cVec create(double* y, const cVec& v) 
  {
    const dcomp dummy=0;
    blablabla::allocator<dcomp> alloc(reinterpret_cast<dcomp*>(y));
    std::vector<dcomp,blablabla::allocator<dcomp> > x(v.size(),dummy,alloc);     
    return x;
  }
  static const cVec create(const double* y, const cVec& v) {return create(const_cast<double*>(y),v);}
};



#endif // _VECTOR_TRAITS_H
