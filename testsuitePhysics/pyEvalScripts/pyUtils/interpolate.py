# -*- coding: utf-8 -*- 

import scipy.interpolate

def interpolate(timeArray,array) :
    return scipy.interpolate.interp1d(timeArray,array) # scipy.interpolate.interp1d(timeArray,array,kind='cubic') takes immense amount of time

def interpolateEVC(evc,i) :
    return interpolate(evc.time,evc[i]) 
