#!/bin/env python2

from cpypyqed import *
import sys

p = parameters.ParameterTable()
pm=mode.ParsPumpedLossy(p)
pq=qbit.ParsPumpedLossy(p)
pe=ParsEvolution(p)
pjc=jaynescummings.Pars(p)

parameters.update(p,sys.argv,'--')

m=mode.make(pm,QMP.IP)
q=qbit.make(pq,QMP.IP)
jc=jaynescummings.make(q,m,pjc)

psi = (qbit.init(pq)*mode.init(pm))
psi.renorm()

evolve(psi, binary.make(jc), pe)
