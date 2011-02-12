# 1particle1mode --configuration 2 --fin 8 --etaeff 0 --eta \(300,-200\) --N 10 --T .1 --Unot 0 --DC -1000 --kappa 100 --pinit \(-0.5,10,0.1\) --dc 1 > free.d

set xra [0:.1]

set term x11 0
set title "IP"
plot "free.dIP" u 1:3 t "photon number", "free.dIP" u 1:4 t "pn dispersion", "free.dIP" u 1:($5**2+$6**2) w l t "|<a>|^2"
set term x11 1
set title "UIP"
plot "free.dUIP" u 1:3 t "photon number", "free.dUIP" u 1:4 t "pn dispersion", "free.dUIP" u 1:($5**2+$6**2) w l t "|<a>|^2"
set term x11 2
set title "Sch"
plot "free.dSch" u 1:3 t "photon number", "free.dSch" u 1:4 t "pn dispersion", "free.dSch" u 1:($5**2+$6**2) w l t "|<a>|^2"
pause -1

etaperzre=-0.168317
etaperzim=-0.316832
alphanotre=-0.5
alphanotim=1

DC=-1000
kappa=100

plot "free.dIP" u 1:(etaperzre+exp(-kappa*$1)*(cos(DC*$1)*(alphanotre-etaperzre)-sin(DC*$1)*(alphanotim-etaperzim))) w l t "exact", "free.dSch" u 1:5 t "<a>.re Sch", "free.dUIP" u 1:5 t "<a>.re UIP", "free.dIP" u 1:5 t "<a>.re IP", "free.dMaster" u 1:5 t "<a>.re Master"
pause -1

plot "free.dIP" u 1:(etaperzim+exp(-kappa*$1)*(cos(DC*$1)*(alphanotim-etaperzim)+sin(DC*$1)*(alphanotre-etaperzre))) w l t "exact", "free.dSch" u 1:6 t "<a>.im Sch", "free.dUIP" u 1:6 t "<a>.im UIP", "free.dIP" u 1:6 t "<a>.im IP", "free.dMaster" u 1:6 t "<a>.im Master"
pause -1

set xra [0:.3]

xnot=-1.57
pnot=10

set ytics ('-' -pi/2, '0' 0, '+' pi/2)
set grid ytics

plot "free.dIP" u 1:(xnot+2*pnot*$1) w l t "exact", "free.dIP" u 1:9 t "<x>", "free.dUIP" u 1:9 t "<x> UIP", "free.dSch" u 1:9 t "<x> Sch", "free.dMaster" u 1:9 t "<x> Master"
pause -1

set ytics ('0' 0, '+' pi/3**.5)

Deltaxnot=0.01
Deltapnot=25

set yra [0:]

plot "free.dIP" u 1:((Deltaxnot+4*Deltapnot*$1**2)**.5) w l t "exact", "free.dIP" u 1:10 t "Delta x", "free.dUIP" u 1:10 t "Delta x UIP", "free.dSch" u 1:10 t "Delta x Sch", "free.dMaster" u 1:10 t "Delta x Master"
pause -1

