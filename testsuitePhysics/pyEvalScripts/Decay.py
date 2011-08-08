#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import regression
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories, argv
from scipy import exp

trajectories=loadTrajectories([

"DecayEn.d"  ,
"DecayMa.d"  ,

"DecayFockEn.d",
"DecayFockMa.d",

],
argv
)


timeExtrema=[0,2]

kappa=2.
gamma=1.

epsIndeces=[0,1,2,1]

def compare(a,l,i,eps) :
    f=lambda t : a*exp(-l*t)
    res=True
    for trajectory, epsi in zip(trajectories,epsIndeces) :
        res&=regression(f,interpolateEVC(trajectory,i),timeExtrema,argv,eps[epsi])
    return res

def comparev(a,l,trajs,te=timeExtrema) :
    f=lambda t : a*exp(-l*t)
    res=True
    for t in trajs :
        res&=regression(f,interpolateEVC(trajectories[t[0]],t[1]),te,argv,t[2])
    return res


eps0=1.5e-8
eps1=7e-8
eps2=1e-7

print compare(.25,2.*gamma,3, [5e-6, 1.4e-8, 3e-6] ) & compare(.346,gamma,4, [3e-5, 1.5e-7, 1.5e-5] ) & compare(.26,gamma,5, [1.5e-5, 1e-8, 7e-6] ) & comparev(1.25,2.*kappa, [ [0,6,eps0], [1,6,eps0], [0,7,eps0], [1,7,eps0] ] ) & comparev(1,kappa, [ [0,8,eps1], [1,8,eps1] ] ) & comparev(-.5,kappa, [ [0,9,eps2], [1,9,eps2] ] ) & comparev(8,2.*kappa, [ [2,6,0.0015], [3,6,3e-9], [2,7,0.005], [3,7,0.0015] ] , [.5,2] )
