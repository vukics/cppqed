etaeff=10000
omega=2*etaeff**.5
xampl=pi*0.01
pampl=etaeff**.5*xampl

set grid ytics

set term x11 0; set ytics ("-" -pampl, "0" 0, "+" pampl)
plot [0:.063][-1.1*pampl:1.1*pampl]  pampl*sin(omega*x), "oscillation10000.d" u 1:3, "oscillation10000_Sch.d" u 1:3 ev 100 w l
set term x11 1; set ytics ("-" -xampl, "0" 0, "+" xampl)
plot [0:.063][-1.1*xampl:1.1*xampl] -xampl*cos(omega*x), "oscillation10000.d" u 1:5, "oscillation10000_Sch.d" u 1:5 ev 100 w l
pause -1

etaeff=1000
omega=2*etaeff**.5
xampl=pi*0.01
pampl=etaeff**.5*xampl


set term x11 0; set ytics ("-" -pampl, "0" 0, "+" pampl)
plot [0:.4][-1.1*pampl:1.1*pampl]  pampl*sin(omega*x), "oscillation1000.d" u 1:3
set term x11 1; set ytics ("-" -xampl, "0" 0, "+" xampl)
plot [0:.4][-1.1*xampl:1.1*xampl] -xampl*cos(omega*x), "oscillation1000.d" u 1:5
pause -1

etaeff=100
omega=2*etaeff**.5
xampl=pi*0.01
pampl=etaeff**.5*xampl


set term x11 0; set ytics ("-" -pampl, "0" 0, "+" pampl)
plot [0:.66][-1.1*pampl:1.1*pampl]  pampl*sin(omega*x), "oscillation100.d" u 1:3
set term x11 1; set ytics ("-" -xampl, "0" 0, "+" xampl)
plot [0:.66][-1.1*xampl:1.1*xampl] -xampl*cos(omega*x), "oscillation100.d" u 1:5
pause -1


#etaeff=10
#omega=2*etaeff**.5
#xampl=pi*0.01
#pampl=etaeff**.5*xampl

#plot "oscillation10.d" u 1:3,  pampl*sin(omega*x); pause -1
#plot "oscillation10.d" u 1:5, -xampl*cos(omega*x); pause -1
