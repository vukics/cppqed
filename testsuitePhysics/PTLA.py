# -*- coding: utf-8 -*-
import pycppqed as qed

from pyUtils.regression import regression

from pyUtils.interpolate import *

PTLA_Ev=qed.load_cppqed("data/PTLA_Ev.d")[0]

PTLA_Ma=qed.load_cppqed("data/QMJ_En.d")[0]

print regression(interpolate(PTLA_Ev.time,(1+PTLA_Ev[2])/2),interpolateEVC(PTLA_Ma,2),PTLA_Ev.time)

