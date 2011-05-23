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


def compare(i,eps) :
    averageArray=trajectories[0][i]
    for trajectory in trajectories[1:6] : averageArray+=trajectory[i]
    for trajectory in trajectories[6:]  : averageArray+=trajectory[i+4][:101]
    averageArray/=9.
    res=True
    for trajectory in trajectories[0:6] :
        reg=regression(interpolate(timeArray,averageArray),interpolateEVC(trajectory,i  ),timeArray)
        # print reg
        res&=reg<eps
    for trajectory in trajectories[6:]  :
        reg=regression(interpolate(timeArray,averageArray),interpolateEVC(trajectory,i+4),timeArray)
        # print reg
        res&=reg<eps
    return res


print compare(2,1e-34) & compare(3,1e-34) & compare(4,1e-34) & compare(5,1e-34)
