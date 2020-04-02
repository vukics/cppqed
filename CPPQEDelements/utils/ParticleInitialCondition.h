// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_UTILS_PARTICLEINITIALCONDITION_H_INCLUDED
#define   CPPQEDELEMENTS_UTILS_PARTICLEINITIALCONDITION_H_INCLUDED

#include "ContainerIO.h"

#include <tuple>


namespace particle {

class InitialCondition
{
public:
  typedef std::tuple<double,double,double,bool> Impl;

  InitialCondition(double x0, double k0, double sig, bool isItInK=false) : impl_(x0,k0,sig,isItInK) {}

  double getX0 () const {return std::get<0>(impl_);}
  double getK0 () const {return std::get<1>(impl_);}
  double getSig() const {return std::get<2>(impl_);}
  bool   isInK () const {return std::get<3>(impl_);}

  const Impl& get() const {return impl_;}
  
  Impl& get() {return const_cast<Impl&>(static_cast<const InitialCondition&>(*this).get());}

private:
  const Impl impl_;
  
};


inline std::ostream& operator<<(std::ostream& os, const InitialCondition& ic) {return os<<ic.get();}
inline std::istream& operator>>(std::istream& is,       InitialCondition& ic) {return is>>ic.get();}



} // particle


#endif // CPPQEDELEMENTS_UTILS_PARTICLEINITIALCONDITION_H_INCLUDED
