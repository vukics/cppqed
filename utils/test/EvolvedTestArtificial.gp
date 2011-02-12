g=1
o=10

plot [0:4] exp(exp(-g*x)*cos(o*x))*cos(exp(-g*x)*sin(o*x)), "EvolvedTestArtificial.dat" u 1:3 "%lf %lf (%lf,%lf) (%lf,%lf)"

set term x11 1

plot [0:4] exp(exp(-g*x)*cos(o*x))*sin(exp(-g*x)*sin(o*x)), "EvolvedTestArtificial.dat" u 1:4 "%lf %lf (%lf,%lf) (%lf,%lf)"

pause -1
