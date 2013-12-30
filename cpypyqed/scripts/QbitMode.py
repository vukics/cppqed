#!/bin/env python2

import sys
if '--debug' in sys.argv:
    from cpypyqed_d import *
    sys.argv.remove('--debug')
else:
    from cpypyqed import *

p = parameters.ParameterTable()
pm=mode.ParsPumpedLossy(p)
pq=qbit.ParsPumpedLossy(p)
pe=ParsEvolution(p)
pjc=jaynescummings.Pars(p)

parameters.update(p,sys.argv,'--')

qmp = QMP.UIP if pe.evol == EM.MASTER or pe.evol == EM.MASTER_FAST else QMP.IP
m=mode.make(pm,qmp)
q=qbit.make(pq,qmp)
jc=jaynescummings.make(q,m,pjc)

psi = (qbit.init(pq)*mode.init(pm))
psi.renorm()

evolve(psi, binary.make(jc), pe)
