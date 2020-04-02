#!/usr/bin/env python3

# the argparse module is only available for Python >= 2.7, this is a version with optparse instead

import sys
from cpypyqed import *
from cpypyqed.tools.passthrough import PassThroughOptionParser

parser = PassThroughOptionParser(add_help_option=False)
parser.add_option('--1p1mconf', dest="conf", help="System configuration code for 1particle1mode",
                    type="int",default=1)
(options,remaining)=parser.parse_args(sys.argv)
conf = options.conf

p  = parameters.ParameterTable()
pe = evolution.Pars(p)

pplm = mode.ParsPumpedLossy(p)
ppp  = particle.ParsPumped(p)
ppci = particlecavity.ParsAlong(p)

qmp = updateWithPicture(p,remaining,'--')

if not conf==1: pplm.delta-=ppci.uNot/(1. if isComplex(ppci.modeCav) else 2.)

myMode     = mode.make(pplm,qmp)
myParticle = particle.make(ppp,qmp)
myPumpedParticle = particle.makePumped(ppp,qmp)

if conf==1:
  particleCavityBase = ParticleOrthogonalToCavity(myMode,myParticle,ppci)
elif conf==2:
  if not abs(pplm.eta) and not ppp.vClass:
    sys.exit("No driving in the system!")
  if not ppp.init.getSig():
    sys.exit("No initial spread specified!")
  particleCavityBase = ParticleAlongCavity(myMode,myParticle,ppci,ppp.vClass)
elif conf==3:
  particleCavityBase = ParticleAlongCavity(myMode,myPumpedParticle,ppci)
elif conf==4:
  particleCavityBase = ParticleAlongCavity(myMode,myPumpedParticle,ppci,0)
else:
  sys.exit("Configuration not recognized!")

psi = (mode.init(pplm)**particle.init(ppp)).normalize()

evolve(psi, binary.make(particleCavityBase), pe)

