#!/bin/env python2

from cpypyqed import *
import sys

p = parameters.ParameterTable()
pm=mode.ParsPumpedLossy(p)
pq=qbit.ParsPumpedLossy(p)
pe=ParsEvolution(p)
pjc=jaynescummings.Pars(p)

parameters.update(p,sys.argv,'--')

qmp = QMP.UIP if pe.evol == EM_MASTER or pe.evol == EM_MASTER_FAST else QMP.IP
m=mode.make(pm,qmp)
q=qbit.make(pq,qmp)
jc=jaynescummings.make(q,m,pjc)

psi = (qbit.init(pq)*mode.init(pm))
psi.renorm()

evolve(psi, binary.make(jc), pe)
