set pointsize .7
set term x11 0
set ytics nomirror
set y2tics
pl [:.2] "QMJ_Int_Ev.d" u 1:5 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 2, "QMJ_Int_Matrix.d" u 1:5 w lp lw 2, "QMJ_Int_Si.d" u 1:5 w lp lw 2, "QMJ_Int_SiSch.d" u 1:5 w lp lw 2, "QMJ_Int_SiUIP.d" u 1:5 w lp lw 2, "QMJ_Int_En.d" u 1:5 w lp lw 2, "QMJ_Int_Ma.d" u 1:5 w lp lw 2, "QMJ_Int_MaSch.d" u 1:5 w lp lw 2, "QMJ_Int_Si.d" u 1:4 w lp lw 2 axes x1y2 
set term x11 1
pl [:.2] "QMJ_Int_Ev.d" u 1:6 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 2, "QMJ_Int_Matrix.d" u 1:6 w lp lw 2, "QMJ_Int_Si.d" u 1:6 w lp lw 2, "QMJ_Int_SiSch.d" u 1:6 w lp lw 2, "QMJ_Int_SiUIP.d" u 1:6 w lp lw 2, "QMJ_Int_En.d" u 1:6 w lp lw 2, "QMJ_Int_Ma.d" u 1:6 w lp lw 2, "QMJ_Int_MaSch.d" u 1:6 w lp lw 2
set term x11 2
pl [:.2] "QMJ_Int_Ev.d" u 1:3 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 2, "QMJ_Int_Matrix.d" u 1:9 w lp lw 2, "QMJ_Int_Si.d" u 1:9 w lp lw 2, "QMJ_Int_SiSch.d" u 1:9 w lp lw 2, "QMJ_Int_SiUIP.d" u 1:9 w lp lw 2, "QMJ_Int_En.d" u 1:9 w lp lw 2, "QMJ_Int_Ma.d" u 1:9 w lp lw 2, "QMJ_Int_MaSch.d" u 1:9 w lp lw 2
set term x11 3
pl [:.2] "QMJ_Int_Ev.d" u 1:4 "%lf %lf (%lf,%lf) (%lf,%lf)" w lp lw 2, "QMJ_Int_Matrix.d" u 1:10 w lp lw 2, "QMJ_Int_Si.d" u 1:10 w lp lw 2, "QMJ_Int_SiSch.d" u 1:10 w lp lw 2, "QMJ_Int_SiUIP.d" u 1:10 w lp lw 2, "QMJ_Int_En.d" u 1:10 w lp lw 2, "QMJ_Int_Ma.d" u 1:10 w lp lw 2, "QMJ_Int_MaSch.d" u 1:10 w lp lw 2
pause -1
