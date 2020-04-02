// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_DO_DISPLAY_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_DO_DISPLAY_TCC_INCLUDED

#include "DO_Display.h"

#include "NegPT.tcc"
#include "Structure.h"

#include "FormDouble.h"


template<int RANK, typename V>
std::ostream&
quantumtrajectory::display_densityoperator::_<RANK,V>::displayKey(std::ostream& os, size_t& i) const 
{
  if (av_) av_->displayKey(os,i); 
  if (negativity_) os<<"# Trajectory\n# "<<i<<". negativity"<<std::endl;
  return os;
}


template<int RANK, typename V>
std::ostream&
quantumtrajectory::display_densityoperator::_<RANK,V>::display(double t, const DensityOperator& rho, std::ostream& os, int precision) const 
{
  structure::display(av_,t,rho,os,precision);
  if (negativity_) os<<'\t'<<FormDouble(precision)(quantumdata::negPT(rho,V()));
  return os;
}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_DO_DISPLAY_TCC_INCLUDED
