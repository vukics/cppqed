#! /usr/bin/env python

import sys
import HCmodule

N=int(sys.argv[1])
hc=HCmodule.hc(N)

print '#include "HermiteCoefficients.H"\nnamespace CppUtils\nnamespace HermiteCoefficients {\n\nsize_t lim=%d;\n' % N

for n in range(N):
    print 'double c%-2d[] = { %12g' % (n,hc[n][0]),
    for k in range(1,n+1): print ', %12g' % hc[n][k],
    print '};'

print '\nstd::vector<double*> Fill() {std::vector<double*> t(%d);' % N,

for n in range(N): print 't[%d]=c%d;' % (n,n),

print 'return t;}\n'

print 'std::vector<double*> c(Fill());\n}\n}'
