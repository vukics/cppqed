// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_TCC_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_TCC_INCLUDED

#include "Pars.tcc"

#include <boost/fusion/sequence/io.hpp>
#include <boost/fusion/include/io.hpp>


namespace multilevel {


template<int NL, typename VP, typename VL>
ParsPumpedLossy<NL,VP,VL>::ParsPumpedLossy(parameters::ParameterTable& p, const std::string& mod)
  : deltas(p.addTitle("PumpedLossyMultiLevel",mod).addMod("deltas",mod,"MultiLevel detunings vector",Levels())),
    etas(p.addMod("etas",mod,"MultiLevel pumps vector",VP())),
    gammas(p.addMod("gammas",mod,"MultiLevel decays vector",VL()))
{
  using namespace boost::fusion;

  deltas=0;
  p.getStream()<<tuple_open('[')<<tuple_close(']');
  std::cout    <<tuple_open('[')<<tuple_close(']');
}


} // multilevel


template<int NL>
std::ostream&
operator<<(std::ostream& os, const blitz::TinyVector<double,NL>& x)
{
  os << '[' << x[0];
  for (int i=1; i < NL; ++i)
    {
      os << ' ' << x[i];
    }
  os << ']';
  return os;
}


template<int NL>
std::istream&
operator>>(std::istream& is,       blitz::TinyVector<double,NL>& x)
{
  using namespace std;

  char sep;
             
  is >> sep;
  BZPRECHECK(sep == '[', "Format error while scanning input TinyVector"
	     << endl << " (expected '[' opening TinyVector)");

  is >> x(0);
  for (int i = 1; i < NL; ++i) {
    BZPRECHECK(!is.bad(), "Premature end of input while scanning TinyVector");
    is >> x(i);
  }
  is >> sep;
  BZPRECHECK(sep == ']', "Format error while scanning input TinyVector"
	     << endl << " (expected ']' closing TinyVector)");
    
  return is;
}


#endif // CPPQEDELEMENTS_FREES_PARSMULTILEVEL_TCC_INCLUDED
