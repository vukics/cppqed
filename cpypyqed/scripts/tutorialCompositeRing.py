#!/usr/bin/env python3
# coding: utf-8
import sys
from cpypyqed import *

p=parameters.ParameterTable()

pe=evolution.Pars(p)

pp=particle.Pars(p)

pmP=mode.ParsLossy(p,"P")
pmM=mode.ParsLossy(p,"M")
ppcP=particlecavity.ParsAlong(p,"P")
ppcM=particlecavity.ParsAlong(p,"M")

ppcP.modeCav=ModeFunctionType.PLUS
ppcM.modeCav=ModeFunctionType.MINUS

qmp=updateWithPicture(p,sys.argv,"--")

part=particle.make(pp ,qmp)
plus =   mode.make(pmP,qmp)
minus=   mode.make(pmM,qmp)

evolve(particle.wavePacket(pp)**mode.init(pmP)**mode.init(pmM),
       makeComposite({(1,0):ParticleAlongCavity(plus ,part,ppcP),
                      (2,0):ParticleAlongCavity(minus,part,ppcM),
                      (1,2,0):ParticleTwoModes(plus,minus,part,ppcP,ppcM)}),
       pe)
