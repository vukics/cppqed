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

p = parameters.ParameterTable()
pe = evolution.Pars(p)
pp = particle.ParsPumped(p)
pm = mode.ParsPumpedLossy(p)
ppc = particlecavity.ParsOrthogonal(p)

parameters.update(p,remaining,'--')

qmp = QMP.UIP if pe.evol == evolution.Method.MASTER or pe.evol == evolution.Method.MASTER_FAST else QMP.IP
myMode     = mode.make(pm,qmp)
myParticle = particle.makePumped(pp,qmp)

poc = ParticleOrthogonalToCavity(myMode,myParticle,ppc)

psi = (mode.init(pm)**particle.init(pp)**particle.init(pp)).normalize()

c = makeComposite({(0,2):poc,(0,1):poc})

# Note that we can delete all the objects which are not needed for evolve.
# Python takes care that the references stored in c and pe are not actually
# deleted. For this to work we have to make sure with_custodian_and_ward is
# used in the wrapper code whenever a reference is held by the returned
# object.
del(p)
del(pp)
del(pm)
del(ppc)
del(myMode)
del(myParticle)
del(poc)
evolve(psi, c, pe)
