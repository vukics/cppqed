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

print helper(2,3e-1) & helper(3,28.) & helper(4,4e-3) & helper(5,7e-3) & helper(6,2e-1) & helper(7,2e-1) & helper(8,6e-4) & helper(9,2e-3) & helper(10,2e-3) & helper(11,2e-3) #& helper(2,1,4,8e-11) & helper(3,0,3,5e-10) & helper(3,1,5,5e-10)

