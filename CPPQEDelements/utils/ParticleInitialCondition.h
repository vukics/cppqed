// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_UTILS_PARTICLEINITIALCONDITION_H_INCLUDED
#define   CPPQEDELEMENTS_UTILS_PARTICLEINITIALCONDITION_H_INCLUDED

#include<boost/tuple/tuple_io.hpp>


namespace particle {

class InitialCondition
{
public:
  typedef boost::tuple<double,double,double,bool> Impl;

  InitialCondition(double x0, double k0, double sig, bool isItInK=false) : impl_(x0,k0,sig,isItInK) {}

  double getX0 () const {return impl_.get<0>();}
  double getK0 () const {return impl_.get<1>();}
  double getSig() const {return impl_.get<2>();}
  bool   isInK () const {return impl_.get<3>();}

  const Impl& get() const {return impl_;}
  
  Impl& get() {return const_cast<Impl&>(static_cast<const InitialCondition&>(*this).get());}

private:
  const Impl impl_;
  
};


inline std::ostream& operator<<(std::ostream& os, const InitialCondition& ic) {return os<<ic.get();}
inline std::istream& operator>>(std::istream& is,       InitialCondition& ic) {return is>>ic.get();}



} // particle


#endif // CPPQEDELEMENTS_UTILS_PARTICLEINITIALCONDITION_H_INCLUDED
