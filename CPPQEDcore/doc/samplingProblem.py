#./PumpedLossyMode_C++QED --cutoff 20 --eta 0 --nTh .1 --T 1 --kappa 10 --dc 1 --dpLimit 0.01 --seed 1002 --evol master > m.d
#./PumpedLossyMode_C++QED --cutoff 20 --eta 0 --nTh .1 --T 1 --kappa 10 --dc 1 --dpLimit 0.01 --evol ensemble --nTraj 10000  > e0.01.d
#./PumpedLossyMode_C++QED --cutoff 20 --eta 0 --nTh .1 --T 1 --kappa 10 --dc 1 --dpLimit 0.1 --evol ensemble --nTraj 10000  > e0.1.d
#./PumpedLossyMode_C++QED --cutoff 20 --eta 0 --nTh .1 --T 1 --kappa 10 --dc 1 --Dt 0.01 --dpLimit 0.1 --evol ensemble --nTraj 10000  > e0.1_dense.d
#./PumpedLossyMode_C++QED --cutoff 20 --eta 0 --nTh .1 --T 2000 --kappa 10 --dc 1 --dpLimit 0.01 --timeAverage --relaxationTime 10 > sNaive.d
#./PumpedLossyMode_C++QED --cutoff 20 --eta 0 --nTh .1 --T 2000 --kappa 10 --dc 0 --Dt 0.1 --dpLimit 0.01 --timeAverage --relaxationTime 10 > s.d

# pl [0:.5] "m.d" u 1:3 w l, "e0.01.d" u 1:3, "e0.1.d" u 1:3, "e0.1_dense.d" u 1:3, 0.103
# pl [104:105] "sNaive.d" u 1:3, "sNaive.d" u 1:2 axes x1y2

import matplotlib
matplotlib.use('SVG')
from matplotlib.pyplot import *
import numpy
rc('text', usetex=True)
rc('font', family='serif')
m=numpy.loadtxt("m.d")
e1=numpy.loadtxt("e0.01.d")
e2=numpy.loadtxt("e0.1.d")
e3=numpy.loadtxt("e0.1_dense.d")
#subplot(2,1,1)
#figure(1, figsize=(6,4))
plot( m[:,0], m[:,2], linewidth=2)
plot( e1[:,0], e1[:,2])
plot( e2[:,0], e2[:,2])
plot( e3[:,0], e3[:,2])
#yscale('log')
xlim(0,.5)
#ylim(5e-5,.05)
xlabel(r'$t$', fontsize=18)
ylabel(r'photon number', fontsize=16)
suptitle(r'The role of dpLimit \& sampling', fontsize=16)
grid()
#subplots_adjust(top=0.8)
#tight_layout()
#show()
savefig('samplingProblem')
