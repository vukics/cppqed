// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "QM_Picture.h"

#include "EvolutionMethod.h"

#include "Pars.h"

#include <iostream>
#include <string>


using namespace std;

ostream& operator<<(ostream& os, QM_Picture qmp)
{
  switch (qmp) {
  case QMP_IP : return os<<"IP" ;
  case QMP_UIP: return os<<"UIP";
  case QMP_SCH:        os<<"Sch";
  }
  return os;
}


istream& operator>>(istream& is, QM_Picture& qmp) 
{
  QM_Picture qmptemp=QMP_SCH;
  string s;

  is>>s;
       if (s=="IP" ) qmptemp=QMP_IP ;
  else if (s=="UIP") qmptemp=QMP_UIP;
  else if (s!="Sch")
    is.clear(ios_base::badbit);

  if (is) qmp=qmptemp;
  return is;
}


QM_Picture& picture::updateWithPicture(parameters::Table& p, int argc, char* argv[], const std::string& prefix)
{
  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);
  parameters::update(p,argc,argv,prefix);
  try {
    const evolution::Method method=dynamic_cast<const parameters::Typed<evolution::Method>&>(p["evol"]).get();
    if (method==evolution::MASTER && qmp==QMP_IP) qmp=QMP_UIP;
  } catch (const parameters::UnrecognisedParameterException&) {}
  return qmp;
}
