#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pyUtils.regression import regressionDiscrete
from pyUtils.loadTrajectories import loadTrajectories, argv

trajectories=loadTrajectories([

"PLM_FT_Ma.d",
"PLM_FT_En.d"

],
argv)


def helper(i,eps) : return regressionDiscrete(trajectories[0],i,trajectories[1],i,argv,eps)

print helper(2,4e-2) & helper(3,6.2) & helper(4,2e-4) & helper(5,4e-4) & helper(6,3e-2) & helper(7,3e-2) & helper(8,4e-4) & helper(9,6e-4) & helper(10,6e-4) & helper(11,5e-4) #& helper(2,1,4,8e-11) & helper(3,0,3,5e-10) & helper(3,1,5,5e-10)

