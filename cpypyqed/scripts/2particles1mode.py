#!/bin/env python2

import sys
if '--debug' in sys.argv:
    from cpypyqed_d import *
    sys.argv.remove('--debug')
else:
    from cpypyqed import *

p = parameters.ParameterTable()
pe = ParsEvolution(p)
pp = particle.ParsPumped(p)
pm = mode.ParsPumpedLossy(p)
ppc = particlecavity.ParsOrthogonal(p)

parameters.update(p,sys.argv,'--')

qmp = QMP.UIP if pe.evol == EM.MASTER or pe.evol == EM.MASTER_FAST else QMP.IP
myMode     = mode.make(pm,qmp)
myParticle = particle.makePumped(pp,qmp)

poc = ParticleOrthogonalToCavity(myMode,myParticle,ppc)

psi = mode.init(pm)*particle.init(pp)*particle.init(pp)

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
