set pointsize .7
set term x11 0
pl [:1] "BiM_Ev.d" u 1:3 "%lf %lf (%lf,%lf)" w lp lw 2, "BiM_Si.d" u 1:5 w lp lw 2
set term x11 1
pl [:1] "BiM_Ev.d" u 1:4 "%lf %lf (%lf,%lf)" w lp lw 2, "BiM_Si.d" u 1:6 w lp lw 2
pause -1
