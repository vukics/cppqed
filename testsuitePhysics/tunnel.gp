set ytics ('-' -pi/2, '0' 0, '+' pi/2)
set grid ytics

set term x11 0

plot [:200][-pi/2:pi/2] "tunnelBH_Wannier.dat" u 1:(-$3) w l lw 2, "tunnel.d" u 1:5 ev 10, "tunnel_Sch.d" u 1:5 ev 10

#set term x11 1

#set ytics auto
#unset grid

#plot [][0:1] "tunnel2particlesBH_Wannier.dat" u 1:($3) w l lw 2, "tunnel2particles.d" u 1:11
pause -1 
