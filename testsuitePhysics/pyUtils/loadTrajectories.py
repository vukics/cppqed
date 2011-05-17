# -*- coding: utf-8 -*-
import pycppqed

def loadTrajectories(files,path) :
    res=[]
    for f in files :
        res.append(pycppqed.load_cppqed(path+"/"+f)[0])
    return res
