# -*- coding: utf-8 -*- 

import scipy.interpolate

def interpolate(timeArray,array) :
    return scipy.interpolate.interp1d(timeArray,array)

def interpolateEVC(evc,i) :
    return interpolate(evc.time,evc[i]) 
