#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob

if len(sys.argv)<4 : 
  print "Arguments:\n1. path to executables\n2. executing script\n3. code for sections"
  sys.exit(1)

print "Executing from", sys.argv[1]

ARGSTraj    = "--dc 0 --Dt .01 --T 2"
ARGSEnsemble= "--nTraj 200 --evol ensemble"

Sections=[
"Two-level atom with different implementations",
"Two-level atom with Heisenberg-Langevin",
"Lossy mode with different implementations",
"Non-interacting composite",
"Pure decay in composite systems",
"Interacting binary",
"Free particle",
"Particle oscillating in classical potential",
"Heavy particle not affected by the cavity field",
"Multilevel system",
"Ring",
"Tunneling",
"NX_CoupledModesElim as oscillator",
"Pure decay in even more composite systems involving alternative decay"
]

for i in range(len(Sections)) : 
  if (2**i & int(sys.argv[3])) != 0 : print Sections[i]