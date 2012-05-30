# -*- coding: utf-8 -*- 

from scipy.integrate import quadrature
from interpolate import interpolate

def regression(f1, f2, timeArray, argv, eps) :
    t0=timeArray[ 0]
    t1=timeArray[-1]
    res=quadrature(lambda t : (f1(t)-f2(t))**2,t0,t1,maxiter=100)[0]
    if len(argv)>2 : print res
    return res<eps


def regressionArrays(a1, a2, timeArray, argv, eps) :
    return regression(interpolate(timeArray,a1),interpolate(timeArray,a2),timeArray,argv,eps)


def regressionDiscrete(a1, i1, a2, i2, argv, eps) :
    r=s=0.
    for v1, v2 in zip(a1[i1],a2[i2]) :
        r+=(v1-v2)**2
        s+=1.
    r/=s
    if len(argv)>2 : print r
    return r<eps


def regressionFunction(a1, i1, f, i2, argv, eps) :
    r=s=0.
    for v1, a2 in zip(a1[i1],a1[i2]) :
        r+=(v1-f(a2))**2
        s+=1.
    r/=s
    if len(argv)>2 : print r
    return r<eps
