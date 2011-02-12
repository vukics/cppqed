#include "QM_Picture.h"

#include<iostream>

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

#include<string>

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
