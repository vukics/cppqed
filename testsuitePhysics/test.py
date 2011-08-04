#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, re
from glob import glob
from os.path import join, split, isfile, isdir
from subprocess import call

bitsPath="bits"

def verifyBitsConsistency() :
  bitsList=glob(join(bitsPath,"[0-9]*"))
  bitsIntegers=[int(split(bit)[1]) for bit in bitsList]
  ret=len(bitsIntegers)
  if ret!=max(bitsIntegers)+1 :
    sys.exit("Inconsistency in bits")
  return ret

numberOfBits=verifyBitsConsistency()
bitsList=[ join(bitsPath,str(i)) for i in range(numberOfBits) ]

sectionTitles=[ re.search("# (.*)\n",open(bit,'r').read()).group(1) for bit in bitsList ]

if len(sys.argv)<5 : 
  print "\nArguments:\n1. code of what to run\n2. path to executables\n3. executing script\n4. data directory\n5. additional parameters (optional)\n\nSections:"
  for i, title in enumerate(sectionTitles) :
    print i, title
  print
  sys.exit(0)

executionCode=sys.argv[1]

print "\nSections to run:\n****************"

necessaryScripts=set()
for i in range(numberOfBits) :
  if (2**i & int(executionCode)) != 0 :
    print sectionTitles[i]
    for line in open(bitsList[i],'r').readlines() :
      m=re.match('l\s*"[^"]*"\s*"([^"\s]*)\s',line)
      if m :
        necessaryScripts.add(m.group(1))
print "\nNecessary scripts:", necessaryScripts, "\n"

scriptPath=sys.argv[2]
for script in necessaryScripts :
  if not os.access(join(scriptPath,script),os.X_OK) :
    sys.exit(script+" does not exist or is not executable!\n")
  
plugin=join("plugins",sys.argv[3]+".sh")
if not isfile(plugin) :
  sys.exit(plugin+" does not exist!\n")

dataDirectory=sys.argv[4]
if not isdir(dataDirectory) :
  sys.exit("No directory "+dataDirectory+"\n")

additionalParameters=" "
if len(sys.argv)>5 :
  additionalParameters+=sys.argv[5]

for i in range(numberOfBits) :
  if (2**i & int(executionCode)) != 0 :
    os.system("bash "+join(bitsPath,"main.sh")+" "+bitsList[i]+" "+scriptPath+" "+plugin+" "+dataDirectory+additionalParameters)

