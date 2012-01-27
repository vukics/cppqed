#set term post enha colo
#set out "temp.eps"

set xra [0:.3]

xnot=-1.57
pnot=10

set ytics ('-' -pi/2, '0' 0, '+' pi/2)
set grid ytics

set term wxt 0

plot "freeParticle.d" u 1:(xnot+2*pnot*$1) w l lw 2 t "exact", "freeParticle.d" u 1:5 w l lw 2 t "<x>", "freeParticle_Sch.d" u 1:5 w l lw 2 t "<x>", "freeParticle1p1m.d" u 1:9 w p lw 2 t "<x>"

set ytics ('0' 0, 'homo' pi/3**.5)

Deltaxnot=0.01
Deltapnot=25

set yra [0:]

set term wxt 1

plot "freeParticle.d" u 1:((Deltaxnot+4*Deltapnot*$1**2)**.5) w l lw 2 t "exact", "freeParticle.d" u 1:6 w l lw 2 t "Delta x", "freeParticle_Sch.d" u 1:6 w l lw 2 t "Delta x", "freeParticle1p1m.d" u 1:10 w p lw 2 t "Delta x"
pause -1

