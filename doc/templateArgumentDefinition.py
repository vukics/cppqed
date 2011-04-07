#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys

file=open(sys.argv[1],'r')

text=file.read()

file.close()

text=text.replace('''Note</p>
<p>Template argument definitions:</p>''','Template argument definitions</p>')

open(sys.argv[1],'w').write(text)
