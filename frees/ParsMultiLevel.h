// -*- C++ -*-
#ifndef   PARS_MULTI_LEVEL_INCLUDED
#define   PARS_MULTI_LEVEL_INCLUDED

#include "MultiLevelFwd.h"

#include "ParsFwd.h"

#include <blitz/tinyvec.h>

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


#include<impl/ParsMultiLevel.tcc>

#endif // PARS_MULTI_LEVEL_INCLUDED
