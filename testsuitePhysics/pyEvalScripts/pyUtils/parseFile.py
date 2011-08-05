# -*- coding: utf-8 -*- 

import numpy

floatingReString=r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)'

complexReString =r'\(\s*'+floatingReString+'\s*,\s*'+floatingReString+'\s*\)'

def parseFile(lineSpec,fileName) :
    return numpy.fromregex(fileName,lineSpec.replace(r'+',r'\s*').replace('f',floatingReString).replace('c',complexReString),numpy.float)
