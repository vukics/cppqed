set pointsize .7
set term x11 0
pl [:1] "PLM_Ev.d" u 1:3 "%lf %lf (%lf,%lf)" w lp lw 2, "PLM_En.d" u 1:5 w lp lw 2, "PLM_Si.d" u 1:5 w lp lw 2, "PLM_SiSch.d" u 1:5 w lp lw 2, "PLM_SiUIP.d" u 1:5 w lp lw 2, "PLM_MaSch.d" u 1:5 w lp lw 2, "PLM_MaUIP.d" u 1:5 w lp lw 2, "QMJ_En.d" u 1:9 w lp lw 2, "QMJ_Ma.d" u 1:9 w lp lw 2, "QMJ_MaSch.d" u 1:9 w lp lw 2
set term x11 1
pl [:1] "PLM_Ev.d" u 1:4 "%lf %lf (%lf,%lf)" w lp lw 2, "PLM_En.d" u 1:6 w lp lw 2, "PLM_Si.d" u 1:6 w lp lw 2, "PLM_SiSch.d" u 1:6 w lp lw 2, "PLM_SiUIP.d" u 1:6 w lp lw 2, "PLM_MaSch.d" u 1:6 w lp lw 2, "PLM_MaUIP.d" u 1:6 w lp lw 2, "QMJ_En.d" u 1:10 w lp lw 2, "QMJ_Ma.d" u 1:10 w lp lw 2, "QMJ_MaSch.d" u 1:10 w lp lw 2
pause -1
