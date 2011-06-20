#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import *
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories
from pyUtils.parseFile import parseFile

trajectories=loadTrajectories([

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

evolvedArray=parseFile("f+f+c","data/PLM_Ev.d")

timeArray=trajectories[0].time


def averageArray(i) :
    res=trajectories[0][i]
    for trajectory in trajectories[1:6] : res+=trajectory[i]
    for trajectory in trajectories[6:]  : res+=trajectory[i+4][:101]
    res/=9.
    return res


def compare(i,eps) :
    res=True
    for trajectory in trajectories[0:6] :
        reg=regression(interpolate(timeArray,averageArray(i)),interpolateEVC(trajectory,i  ),timeArray)
        # print reg
        res&=reg<eps
    for trajectory in trajectories[6:]  :
        reg=regression(interpolate(timeArray,averageArray(i)),interpolateEVC(trajectory,i+4),timeArray)
        # print reg
        res&=reg<eps
    return res


alphaRe=evolvedArray[:,2]
alphaIm=evolvedArray[:,3]

print compare(2,1e-34) & compare(3,1e-34) & compare(4,1e-34) & compare(5,1e-34) & (regressionArrays(averageArray(4),alphaRe,timeArray)<1e-8) & (regressionArrays(averageArray(5),alphaIm,timeArray)<1e-7) & (regressionArrays(averageArray(2),alphaRe*alphaRe+alphaIm*alphaIm,timeArray)<1e-8)
