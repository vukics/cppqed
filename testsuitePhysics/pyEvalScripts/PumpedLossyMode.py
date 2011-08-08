#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import *
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories, argv, decidePath, join
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
argv)

evolvedArray=parseFile("f+f+c",join(decidePath(argv),"PLM_Ev.d"))

timeArray=trajectories[0].time


def averageArray(i) : # only the Master ones
    return (trajectories[4][i]+trajectories[5][i]+trajectories[7][i+4][:101]+trajectories[8][i+4][:101])/4.


def compare(i,eps) :
    res=True
    for trajectory in trajectories[0:6] : res&=regression(interpolate(timeArray,averageArray(i)),interpolateEVC(trajectory,i  ),timeArray,argv,eps)
    for trajectory in trajectories[6:]  : res&=regression(interpolate(timeArray,averageArray(i)),interpolateEVC(trajectory,i+4),timeArray,argv,eps)
    return res


def helper(i,array,eps) : return regressionArrays(averageArray(i),array,timeArray,argv,eps)


alphaRe=evolvedArray[:,2]
alphaIm=evolvedArray[:,3]

print compare(2,1e-34) & compare(3,1e-34) & compare(4,1e-34) & compare(5,1e-34) & helper(4,alphaRe,5e-9) & helper(5,alphaIm,5e-8) & helper(2,alphaRe*alphaRe+alphaIm*alphaIm,9e-9)
