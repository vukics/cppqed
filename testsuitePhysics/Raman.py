#! /usr/bin/env python

class MyError(Exception):
    "Exception class"

import sys
import os

if len(sys.argv)<6 : 
    print "delta20 Re(eta0) Im(eta0) Re(eta1) Im(eta1)"
    sys.exit(0)


delta20=float  (sys.argv[1])
eta0   =complex(float(sys.argv[2]),float(sys.argv[3]))
eta1   =complex(float(sys.argv[4]),float(sys.argv[5]))

etat=eta0.conjugate()*eta1/complex(0,delta20);

print "--deltaA", (abs(eta1)**2-abs(eta0)**2)/delta20, "--etat '("+str(etat.real)+","+str(etat.imag)+")'" 
