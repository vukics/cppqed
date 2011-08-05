# -*- coding: utf-8 -*-
import pycppqed
from os.path import join
from sys import argv

def decidePath(argv) :
    if len(argv)<2 : return "."
    else :           return argv[1]


def loadTrajectories(files,argv) :
    res=[]
    for f in files : res.append(pycppqed.load_cppqed(join(decidePath(argv),f))[0])
    return res
