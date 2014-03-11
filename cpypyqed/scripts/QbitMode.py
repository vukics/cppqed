#!/usr/bin/env python2

import sys
import argparse

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--debug', action='store_true')
(args,remaining)=parser.parse_known_args(sys.argv)

if vars(args)['debug']:
    from cpypyqed_d import *
else:
    from cpypyqed import *

p=parameters.ParameterTable()
pm=mode.ParsPumpedLossy(p)
pq=qbit.ParsPumpedLossy(p)
pe=evolution.Pars(p)
pjc=jaynescummings.Pars(p)

parameters.update(p,remaining,'--')

qmp = QMP.UIP if pe.evol == evolution.Method.MASTER or pe.evol == evolution.Method.MASTER_FAST else QMP.IP
m=mode.make(pm,qmp)
q=qbit.make(pq,qmp)
jc=jaynescummings.make(q,m,pjc)

psi = (qbit.init(pq)**mode.init(pm)).normalize()

evolve(psi, binary.make(jc), pe)

