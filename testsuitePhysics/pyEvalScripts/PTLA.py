#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import regression
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories, argv

trajectories=loadTrajectories([
# The "exact" ones:

"PTLA_Ev.d"  ,
"PTLA_Ma.d"  ,

"PLQ_MaSch.d",
"PLQ_MaUIP.d",

"QMJ_Ma.d"   ,
"QMJ_MaSch.d",

# The less "exact" ones:

"PTLA_En.d"  ,

"PLQ_En.d"   ,
"PLQ_EnSch.d",
"PLQ_EnUIP.d",

# One more:

"QMJ_En.d"
],
argv
)

epsIndeces=[0,0,0,0,0,1,1,1,1,2]

timeArray=trajectories[0].time

def compare(array,i,eps) :
    res=True
    for trajectory, epsi in zip(trajectories[1:],epsIndeces) :
        res&=regression(interpolate(timeArray,array),interpolateEVC(trajectory,i),timeArray,argv,eps[epsi])
    return res

print compare((1+trajectories[0][2])/2,2, [ 1.5e-7, .0008, .0003 ] ) & compare((1-trajectories[0][2])/2,3, [ 1.3e-7, .0008, .0003 ] ) & compare(trajectories[0][3]/2,4, [ 2.2e-10, .0002, 4e-5 ] ) & compare(trajectories[0][4]/2,5, [ 1.2e-7, .0007, .00015 ] )
