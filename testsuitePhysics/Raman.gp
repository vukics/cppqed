set term x11

set pointsize .9

plot [][0:1] "RamanTwoLevel.d" u 1:3 w lp lw 3, "RamanRegression.d" u 1:3 w lp lw 3, "RamanCavity.d" u 1:3 ev 10 w lp lw 3

pause -1

plot [][0:1] "RamanTwoLevel.d" u 1:4 w lp lw 3, "RamanRegression.d" u 1:4 w lp lw 3, "RamanCavity.d" u 1:4 ev 10 w lp lw 3

pause -1
