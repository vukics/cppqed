set xrange [0:0.05]

#plot "confAZ_Sin_K1.d" u 1:10, "confAZ_Sin_K2.d" u 1:10, "confAZ_Cos.d" u 1:10, "confAZ_Minus.d" u 1:10, "confAZ_Plus.d" u 1:10; pause -1
#plot "confAZ_Sin_K1.d" u 1:9, "confAZ_Sin_K2.d" u 1:9, "confAZ_Cos.d" u 1:9, "confAZ_Minus.d" u 1:9, "confAZ_Plus.d" u 1:9; pause -1

Unot=2e7
etaeff=0.1
kappa=1e5

A=(Unot*etaeff)**.5/1e5/2

fre(x)=sin(x)
fim(x)=0

alphare(x)=A*fre(x)-A*fim(x)
alphaim(x)=A*fre(x)+A*fim(x)

set term x11 0
plot [0:0.01] "confAZ_Sin_K1.d" u 1:5, "confAZ_Sin_K1.d" u 1:(alphare($9)) w l lw 2
set term x11 1
plot [0:0.01] "confAZ_Sin_K1.d" u 1:6, "confAZ_Sin_K1.d" u 1:(alphaim($9)) w l lw 2

fre(x)=sin(2*x)
fim(x)=0

set term x11 2
plot [0:0.01] "confAZ_Sin_K2.d" u 1:5, "confAZ_Sin_K2.d" u 1:(alphare($9)) w l lw 2
set term x11 3
plot [0:0.01] "confAZ_Sin_K2.d" u 1:6, "confAZ_Sin_K2.d" u 1:(alphaim($9)) w l lw 2

fre(x)=cos(x)
fim(x)=sin(x)

set term x11 4
plot [0:0.01] "confAZ_Plus.d" u 1:5, "confAZ_Plus.d" u 1:(alphare($9)) w l lw 2
set term x11 5
plot [0:0.01] "confAZ_Plus.d" u 1:6, "confAZ_Plus.d" u 1:(alphaim($9)) w l lw 2

fre(x)=cos(x)
fim(x)=-sin(x)

set term x11 6
plot [0:0.01] "confAZ_Minus.d" u 1:5, "confAZ_Minus.d" u 1:(alphare($9)) w l lw 2
set term x11 7
plot [0:0.01] "confAZ_Minus.d" u 1:6, "confAZ_Minus.d" u 1:(alphaim($9)) w l lw 2
pause -1
