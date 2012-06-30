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

print helper(2,3e-1) & helper(3,3e1) & helper(4,3e-3) & helper(5,5e-3) & helper(6,1.5e-1) & helper(7,1.5e-1) & helper(8,9e-4) & helper(9,1.5e-3) & helper(10,1.5e-3) & helper(11,1e-3) #& helper(2,1,4,8e-11) & helper(3,0,3,5e-10) & helper(3,1,5,5e-10)

