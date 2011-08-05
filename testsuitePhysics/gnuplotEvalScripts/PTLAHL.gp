set pointsize .7
set term x11 0
pl [:.1] "PTLAHL_Ev.d" u 1:3 w lp lw 2, "PTLAHL_Si.d" u 1:5 w lp lw 2
set term x11 1
pl [:.1] "PTLAHL_Ev.d" u 1:4 w lp lw 2, "PTLAHL_Si.d" u 1:6 w lp lw 2
pause -1

#pause -1
