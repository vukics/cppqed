#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import regression
from pyUtils.interpolate import *
from pyUtils.loadTrajectories import loadTrajectories

trajectories=loadTrajectories([

"PTLAHL_Ev.d",
"PTLAHL_Si.d"

],
"data")


timeArray=trajectories[0].time


def averageArray(i) :
    return (trajectories[0][i]+trajectories[1][i+2])/2

print (regression(interpolate(timeArray,averageArray(2)),interpolateEVC(trajectories[0],2),timeArray)<2e-11) & (regression(interpolate(timeArray,averageArray(2)),interpolateEVC(trajectories[1],4),timeArray)<2e-11) & (regression(interpolate(timeArray,averageArray(3)),interpolateEVC(trajectories[0],3),timeArray)<1e-11) & (regression(interpolate(timeArray,averageArray(3)),interpolateEVC(trajectories[1],5),timeArray)<1e-11)

