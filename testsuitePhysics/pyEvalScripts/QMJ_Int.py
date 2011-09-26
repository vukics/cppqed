#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import *
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories, argv, decidePath, join
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
argv)

evolvedArray=parseFile("f+f+c+c",join(decidePath(argv),"QMJ_Int_Ev.d"))

timeArray=trajectories[0].time


def averageArray(i) :
    return (trajectories[5][i]+trajectories[6][i])/2.


def compare(i,eps) :
    res=True
    for trajectory in trajectories :
        res&=regression(interpolate(timeArray,averageArray(i)),interpolateEVC(trajectory,i),timeArray,argv,eps)
    return res

def helper(i,j,eps) : return regressionArrays(averageArray(i),evolvedArray[:,j],timeArray,argv,eps)

print compare(4,1e-40) & compare(5,1e-40) & compare(8,1e-40) & compare(9,1e-40) & helper(4,4,3e-8) & helper(5,5,6e-9) & helper(8,2,7e-6) & helper(9,3,5e-6)
