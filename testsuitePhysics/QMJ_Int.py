#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import *
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories
from pyUtils.parseFile import parseFile

trajectories=loadTrajectories([

"QMJ_Int_Matrix.d",
"QMJ_Int_Si.d",
"QMJ_Int_SiSch.d",
"QMJ_Int_SiUIP.d",
"QMJ_Int_En.d",
"QMJ_Int_Ma.d",
"QMJ_Int_MaSch.d",
"QMJ_Int_Si.d",

],
"data")

evolvedArray=parseFile("f+f+c+c","data/QMJ_Int_Ev.d")

timeArray=trajectories[0].time


def averageArray(i) :
    res=trajectories[0][i]
    for trajectory in trajectories[1:] : res+=trajectory[i]
    res/=8.
    return res


def compare(i,eps) :
    res=True
    for trajectory in trajectories :
        reg=regression(interpolate(timeArray,averageArray(i)),interpolateEVC(trajectory,i),timeArray)
        #print reg
        res&=reg<eps
    return res


print compare(4,1e-37) & compare(5,3e-38) & compare(8,1e-40) & compare(9,1e-40) & (regressionArrays(averageArray(4),evolvedArray[:,4],timeArray)<3e-8) & (regressionArrays(averageArray(5),evolvedArray[:,5],timeArray)<6e-9) & (regressionArrays(averageArray(8),evolvedArray[:,2],timeArray)<7e-6) & (regressionArrays(averageArray(9),evolvedArray[:,3],timeArray)<5e-6)
