#!/usr/bin/env python3

import sys
from cpypyqed import *

p=parameters.ParameterTable()

pe=evolution.Pars(p)

pq=qbit.ParsPumpedLossy(p)
pm=mode.ParsPumpedLossy(p)
pjc=jaynescummings.Pars(p)

qmp=updateWithPicture(p,sys.argv,'--')

evolve(qbit.init(pq)**mode.init(pm),
       binary.make( jaynescummings.make( qbit.make(pq,qmp),
                                         mode.make(pm,qmp),
                                         pjc ) ),
       pe)

