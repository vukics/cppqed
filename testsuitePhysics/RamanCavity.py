#! /usr/bin/env python

class MyError(Exception):
    "Exception class"

import sys
import os

if len(sys.argv)<5 : 
    print "Re(eta1) Im(eta1)"
    sys.exit(0)


eta1=complex(float(sys.argv[3]),float(sys.argv[4]))

alpha=complex(11,-9)

z=complex(-float(sys.argv[2]),float(sys.argv[1]))

g=eta1/alpha

eta=-alpha*z

print "--eta '("+str(eta.real)+","+str(eta.imag)+")' --gs '(("+str(g.real)+","+str(g.imag)+"))'" 
