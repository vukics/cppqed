// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED

#include "MultiLevel_Fwd.h"

#include "ParsFwd.h"

#include "BlitzTiny.h"


namespace multilevel {


template<int NL, typename VP, typename VL>
struct ParsPumpedLossy
{
  typedef blitz::TinyVector<double,NL> Levels;

  Levels& deltas;
  VP& etas;
  VL& gammas;

  ParsPumpedLossy(parameters::ParameterTable&, const std::string& ="");

};


} // multilevel


template<int NL>
std::ostream&
operator<<(std::ostream&, const blitz::TinyVector<double,NL>&);

template<int NL>
std::istream&
operator>>(std::istream&,       blitz::TinyVector<double,NL>&);


#include<ParsMultiLevel.tcc>

#endif // CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED
