set pointsize .9

set term x11 0
pl "Ring_Ev.d" u 1:3 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "Ring.d" u 1:9
set term x11 1
pl "Ring_Ev.d" u 1:4 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "Ring.d" u 1:10
set term x11 2
pl "Ring_Ev.d" u 1:5 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "Ring.d" u 1:13
set term x11 3
pl "Ring_Ev.d" u 1:6 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "Ring.d" u 1:14
pause -1

set term x11 0
pl "RingSC_Ev.d" u 1:3 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "RingSC.d" u 1:9
set term x11 1
pl "RingSC_Ev.d" u 1:4 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "RingSC.d" u 1:10
set term x11 2
pl "RingSC_Ev.d" u 1:5 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "RingSC.d" u 1:13
set term x11 3
pl "RingSC_Ev.d" u 1:6 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "RingSC.d" u 1:14
pause -1

set term x11 0
pl "RingS2S_Ev.d" u 1:3 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "RingS2S.d" u 1:9
set term x11 1
pl "RingS2S_Ev.d" u 1:4 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "RingS2S.d" u 1:10
set term x11 2
pl "RingS2S_Ev.d" u 1:5 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "RingS2S.d" u 1:13
set term x11 3
pl "RingS2S_Ev.d" u 1:6 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 3, "RingS2S.d" u 1:14
pause -1


