set xrange [0:0.05]

#plot "confAZ_Sin_K1.d" u 1:10, "confAZ_Sin_K2.d" u 1:10, "confAZ_Cos.d" u 1:10, "confAZ_Minus.d" u 1:10, "confAZ_Plus.d" u 1:10; pause -1
#plot "confAZ_Sin_K1.d" u 1:9, "confAZ_Sin_K2.d" u 1:9, "confAZ_Cos.d" u 1:9, "confAZ_Minus.d" u 1:9, "confAZ_Plus.d" u 1:9; pause -1

Unot=-6e4
etaeff=-10
kappa=1e5
DC=-1e5

R=Unot/kappa

A=(Unot*etaeff)**.5/kappa

fre(x)=sin(x)
fim(x)=0

aux(x)=1+R*(fre(x)**2+fim(x)**2)
auxx(x)=1+aux(x)**2

auxre(x)=aux(x)/auxx(x)
auxim(x)=1/auxx(x)

alphare(x)=A*(fre(x)*auxre(x)-fim(x)*auxim(x))
alphaim(x)=A*(fre(x)*auxim(x)+fim(x)*auxre(x))


set term x11 0
plot [0:0.01] "confAX_Sin_K1.d" u 1:5, "confAX_Sin_K1.d" u 1:(alphare($9)) w l lw 2
set term x11 1
plot [0:0.01] "confAX_Sin_K1.d" u 1:6, "confAX_Sin_K1.d" u 1:(alphaim($9)) w l lw 2

fre(x)=sin(2*x)
fim(x)=0

set term x11 2
plot [0:0.01] "confAX_Sin_K2.d" u 1:5, "confAX_Sin_K2.d" u 1:(alphare($9)) w l lw 2
set term x11 3
plot [0:0.01] "confAX_Sin_K2.d" u 1:6, "confAX_Sin_K2.d" u 1:(alphaim($9)) w l lw 2

fre(x)=cos(x)
fim(x)=0

set term x11 4
plot [0:0.01] "confAX_Cos.d" u 1:5, "confAX_Cos.d" u 1:(alphare($9)) w l lw 2
set term x11 5
plot [0:0.01] "confAX_Cos.d" u 1:6, "confAX_Cos.d" u 1:(alphaim($9)) w l lw 2

fre(x)=cos(x)
fim(x)=-sin(x)

set term x11 6
plot [0:0.01] "confAX_Minus.d" u 1:5, "confAX_Minus.d" u 1:(alphare($9)) w l lw 2
set term x11 7
plot [0:0.01] "confAX_Minus.d" u 1:6, "confAX_Minus.d" u 1:(alphaim($9)) w l lw 2

fre(x)=cos(x)
fim(x)=sin(x)

set term x11 8
plot [0:0.01] "confAX_Plus.d" u 1:5, "confAX_Plus.d" u 1:(alphare($9)) w l lw 2
set term x11 9
plot [0:0.01] "confAX_Plus.d" u 1:6, "confAX_Plus.d" u 1:(alphaim($9)) w l lw 2
pause -1
