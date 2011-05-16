#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob

if len(sys.argv)<4 : 
  print "Arguments:\n1. path to executables\n2. executing script\n3. code for sections"
  sys.exit(1)

path=sys.argv[1]+'/'

print "Executing from", path, "\n"

ARGSTraj    ="--dc 0 --Dt .01 --T 2 "
ARGSEnsemble="--nTraj 200 --evol ensemble "

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


def run(dataFileName,commandLine) :
  print dataFileName
  os.system(path+commandLine+' >'+dataFileName)

def runSection0() :

  ARGS0Base="--etat '(1,3)' --gamma 1.2 --deltaA -1 --qbitInit '(.4,.3)' "

  ARGS0=ARGSTraj+ARGS0Base

  run("PTLA_Ev.d  ","PTLA_Evolved    "+ARGS0)
  run("PTLA_Ma.d  ","PTLA_C++QED     "+ARGS0+"              --evol master")

  run("PTLA_En.d  ","PTLA_C++QED     "+ARGS0+"              $ARGSEnsemble")

  run("PLQ_En.d   ","PumpedLossyQbit "+ARGS0+"              $ARGSEnsemble")
  run("PLQ_EnSch.d","PumpedLossyQbit "+ARGS0+"--picture Sch $ARGSEnsemble")
  run("PLQ_EnUIP.d","PumpedLossyQbit "+ARGS0+"--picture UIP $ARGSEnsemble")
  run("PLQ_MaSch.d","PumpedLossyQbit "+ARGS0+"--picture Sch --evol master")
  run("PLQ_MaUIP.d","PumpedLossyQbit "+ARGS0+"--picture UIP --evol master")


SectionFunctions=[
  runSection0,
]


for i in range(len(Sections)) : 
  if (2**i & int(sys.argv[3])) != 0 :
    print Sections[i]
    SectionFunctions[i]()





