#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import regression
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories

trajectories=loadTrajectories([

#"PLM_Ev.d", # u 1:3 "%lf %lf (%lf,%lf)" w lp lw 2, 
"PLM_En.d",
"PLM_Si.d",
"PLM_SiSch.d",
"PLM_SiUIP.d",
"PLM_MaSch.d",
"PLM_MaUIP.d",

"QMJ_En.d",
"QMJ_Ma.d",
"QMJ_MaSch.d"

],
"data")

timeArray=trajectories[0].time

averageArray=trajectories[0][4]
for trajectory in trajectories[1:5] : averageArray+=trajectory[4]
for trajectory in trajectories[6:]  : averageArray+=trajectory[8][:101]
averageArray/=9

def compare(i) :
    # res=True
    for trajectory in trajectories :
        #res&=
        print regression(interpolate(timeArray,averageArray),interpolateEVC(trajectory,i),timeArray)
    # return res

compare(4)

# print compare((1+trajectories[0][2])/2,2, [ 1.5e-7, .0008, .0003 ] ) & compare((1-trajectories[0][2])/2,3, [ 1.5e-7, .0008, .0003 ] ) & compare(trajectories[0][3]/2,4, [ 3e-10, .0002, 4e-5 ] ) & compare(trajectories[0][4]/2,5, [ 1.2e-7, .0007, .00015 ] )
