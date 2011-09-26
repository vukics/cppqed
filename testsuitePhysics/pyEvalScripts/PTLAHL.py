#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import regression
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories, argv

trajectories=loadTrajectories([

"PTLAHL_Ev.d",
"PTLAHL_Si.d"

],
argv)


timeArray=trajectories[0].time


def averageArray(i) : return (trajectories[0][i]+trajectories[1][i+2])/2

def helper(i,j,k,eps) : return regression(interpolate(timeArray,averageArray(i)),interpolateEVC(trajectories[j],k),timeArray,argv,eps)

print helper(2,0,2,8e-11) & helper(2,1,4,8e-11) & helper(3,0,3,5e-10) & helper(3,1,5,5e-10)

