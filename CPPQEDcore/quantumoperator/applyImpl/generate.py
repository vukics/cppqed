#! /usr/bin/env python
# -*- coding: utf-8 -*-
from sys import argv

rank=int(argv[1])

j=0

def inc(j) :
    j+=1


def recurse(j, i, s0, s1, s2) :
    if (i<rank) :
        s=",I"+str(i)
        return recurse(recurse(recurse(j, i+1, s0+s+" ", s1+s+" ", s2+s+" "), i+1, s0+s+"h", s1+s+"l", s2+s+"l"), i+1, s0+s+"l", s1+s+"l", s2+s+"h")
    else :
        # print j, s0, s1, s2
        print "  if (a("+str(j)+").size()) dpsidt("+s0+")+=a("+str(j)+")("+s1+")*psi("+s2+");"
        return j+1

print """#include "Tridiagonal.h"


namespace quantumoperator {


template<>
void Tridiagonal<"""+argv[1]+""">::apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  using blitz::Range;
  
  const Diagonals & a=diagonals_  ;
  const Dimensions& K=differences_;

  Range"""

for i in range(rank) :
    s=str(i)
    print """    I"""+s+" (0,psi.ubound("+s+""")),
    I"""+s+"l(0,psi.ubound("+s+")-int(K("+s+"""))),
    I"""+s+"h(I"+s+"l+int(K("+s+")))"+("," if i<rank-1 else ";\n")

recurse(recurse(recurse(0,1,"I0 ","I0 ","I0 "),1,"I0h","I0l","I0l"),1,"I0l","I0l","I0h")

print """
}


} // quantumoperator
"""
