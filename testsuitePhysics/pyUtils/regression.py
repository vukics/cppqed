# -*- coding: utf-8 -*- 

import scipy.integrate

def regression(f1, f2, timeArray) :
    return scipy.integrate.quadrature(lambda t : ((f1(t)-f2(t))/(f1(t)+f2(t)))**2,timeArray[0],timeArray[-1],maxiter=100)[0]/(timeArray[-1]-timeArray[0])

