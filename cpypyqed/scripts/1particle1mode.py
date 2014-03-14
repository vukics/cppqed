#!/usr/bin/env python2

import sys
import argparse
from cpypyqed import *

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--1p1mconf', help="System configuration code for 1particle1mode",
                    type=int,default=1)
(args,remaining)=parser.parse_known_args(sys.argv)
conf = vars(args)['1p1mconf']

p  = parameters.ParameterTable()
pe = evolution.Pars(p)

pplm = mode.ParsPumpedLossy(p)
ppp  = particle.ParsPumped(p)
ppci = particlecavity.ParsAlong(p)

parameters.update(p,remaining,'--')

qmp = QMP.UIP if pe.evol == evolution.Method.MASTER or pe.evol == evolution.Method.MASTER_FAST else QMP.IP

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

