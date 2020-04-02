#!/usr/bin/env python3

import sys
from cpypyqed import *

p=parameters.ParameterTable()

pe=evolution.Pars(p)
pm=mode.ParsPumpedLossy(p)

pe.evol=evolution.Method.MASTER
pm.cutoff=30

parameters.update(p,sys.argv,'--')

evolve(mode.init(pm),mode.make(pm,QMP.UIP),pe)
