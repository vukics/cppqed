#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import regressionFunction
from pyUtils.loadTrajectories import loadTrajectories, argv

trajectories=loadTrajectories([

"freeParticle.d",
"freeParticle_Sch.d"

],
argv)

xnot=-1.57
pnot=10

Deltaxnot=0.01
Deltapnot=25

def xC(t) : return xnot+2*pnot*t

def varXC(t) : return (Deltaxnot+4*Deltapnot*t**2)**.5

eps1=5e-5
eps2=1.5e-6

print regressionFunction(trajectories[0][:,0:21],4,xC,0,argv,eps1) & regressionFunction(trajectories[1][:,0:21],4,xC,0,argv,eps1) & regressionFunction(trajectories[0][:,0:21],5,varXC,0,argv,eps2) & regressionFunction(trajectories[1][:,0:21],5,varXC,0,argv,eps2)
