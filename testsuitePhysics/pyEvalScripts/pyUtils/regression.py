# -*- coding: utf-8 -*- 

from scipy.integrate import quadrature
from interpolate import interpolate

def regression(f1, f2, timeArray, argv) :
    t0=timeArray[ 0]
    t1=timeArray[-1]
    res=quadrature(lambda t : (f1(t)-f2(t))**2,t0,t1,maxiter=100)[0]#/quadrature(lambda t : f1(t)*f2(t),t0,t1,maxiter=100)[0]
    # Normalization does not seem to help at all.
    if len(argv)>2 : print res
    return res


def regressionArrays(a1, a2, timeArray, argv) :
    return regression(interpolate(timeArray,a1),interpolate(timeArray,a2),timeArray,argv)


