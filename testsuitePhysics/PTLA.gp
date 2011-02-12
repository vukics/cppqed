set pointsize .7
set term x11 0
pl [:2] "PTLA_Ev.d" u 1:(1+$3)/2 w lp lw 2, "PTLA_Ma.d" u 1:3 w lp lw 2, "PTLA_En.d" u 1:3 w lp lw 2, "PLQ_En.d" u 1:3 w lp lw 2, "PLQ_EnSch.d" u 1:3 w lp lw 2, "PLQ_EnUIP.d" u 1:3 w lp lw 2, "PLQ_MaSch.d" u 1:3 w lp lw 2, "PLQ_MaUIP.d" u 1:3 w lp lw 2, "QMJ_En.d" u 1:3 w lp lw 2, "QMJ_Ma.d" u 1:3 w lp lw 2, "QMJ_MaSch.d" u 1:3 w lp lw 2
set term x11 1
pl [:2] "PTLA_Ev.d" u 1:(1-$3)/2 w lp lw 2, "PTLA_Ma.d" u 1:4 w lp lw 2, "PTLA_En.d" u 1:4 w lp lw 2, "PLQ_En.d" u 1:4 w lp lw 2, "PLQ_EnSch.d" u 1:4 w lp lw 2, "PLQ_EnUIP.d" u 1:4 w lp lw 2, "PLQ_MaSch.d" u 1:4 w lp lw 2, "PLQ_MaUIP.d" u 1:4 w lp lw 2, "QMJ_En.d" u 1:4 w lp lw 2, "QMJ_Ma.d" u 1:4 w lp lw 2, "QMJ_MaSch.d" u 1:4 w lp lw 2
set term x11 2
pl [:2] "PTLA_Ev.d" u 1:($4)/2 w lp lw 2, "PTLA_Ma.d" u 1:5 w lp lw 2, "PTLA_En.d" u 1:5 w lp lw 2, "PLQ_En.d" u 1:5 w lp lw 2, "PLQ_EnSch.d" u 1:5 w lp lw 2, "PLQ_EnUIP.d" u 1:5 w lp lw 2, "PLQ_MaSch.d" u 1:5 w lp lw 2, "PLQ_MaUIP.d" u 1:5 w lp lw 2, "QMJ_En.d" u 1:5 w lp lw 2, "QMJ_Ma.d" u 1:5 w lp lw 2, "QMJ_MaSch.d" u 1:5 w lp lw 2
set term x11 3
pl [:2] "PTLA_Ev.d" u 1:($5)/2 w lp lw 2, "PTLA_Ma.d" u 1:6 w lp lw 2, "PTLA_En.d" u 1:6 w lp lw 2, "PLQ_En.d" u 1:6 w lp lw 2, "PLQ_EnSch.d" u 1:6 w lp lw 2, "PLQ_EnUIP.d" u 1:6 w lp lw 2, "PLQ_MaSch.d" u 1:6 w lp lw 2, "PLQ_MaUIP.d" u 1:6 w lp lw 2, "QMJ_En.d" u 1:6 w lp lw 2, "QMJ_Ma.d" u 1:6 w lp lw 2, "QMJ_MaSch.d" u 1:6 w lp lw 2
pause -1
#pl "PTLA_Ev.d" u 1:2 w lp lw 2, "PTLA_Ma.d" u 1:2 w lp lw 2, "PTLA_En.d" u 1:2 w lp lw 2, "PLQ_En.d" u 1:2 w lp lw 2, "PLQ_EnSch.d" u 1:2 w lp lw 2, "PLQ_EnUIP.d" u 1:2 w lp lw 2, "PLQ_MaSch.d" u 1:2 w lp lw 2, "PLQ_MaUIP.d" u 1:2 w lp lw 2
#pause -1
