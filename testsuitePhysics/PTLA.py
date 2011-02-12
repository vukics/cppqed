# -*- coding: utf-8 -*-
import pycppqed as qed

exact=qed.load_cppqed("PTLA_Ev.d")
toBeVerified=qed.load_cppqed("QMJ_En.d")

evs=toBeVerified[0]

import scipy.interpolate, scipy.integrate

exactInterp=scipy.interpolate.interp1d(exact[0].time,exact[0][3,:]/2.)

toBeVerifiedInterp=scipy.interpolate.interp1d(evs.time,evs[4,:])

print scipy.integrate.quadrature(lambda t : (exactInterp(t)-toBeVerifiedInterp(t))**2,exact[0].time[0],exact[0].time[-1])/(exact[0].time[-1]-exact[0].time[0])

print ((exact[0][3,:]/2.-evs[4,:])**2).mean()
