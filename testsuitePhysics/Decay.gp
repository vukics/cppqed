set pointsize .7

d(i,k,x)=i*exp(-k*x)

kappa=2
gamma=1

e1(i,x)=d(i,  kappa,x)
e2(i,x)=d(i,2*kappa,x)

f1(i,x)=d(i,  gamma,x)
f2(i,x)=d(i,2*gamma,x)

set term x11 0; cc=4; pl [:] "DecayEn.d" u 1:cc w lp lw 2, "DecayMa.d" u 1:cc w lp lw 2, "DecayFockEn.d" u 1:cc w lp lw 2, "DecayFockMa.d" u 1:cc w lp lw 2, f2(.25,x)
set term x11 1; cc=5; pl [:] "DecayEn.d" u 1:cc w lp lw 2, "DecayMa.d" u 1:cc w lp lw 2, "DecayFockEn.d" u 1:cc w lp lw 2, "DecayFockMa.d" u 1:cc w lp lw 2, f1(.346,x)
set term x11 2; cc=6; pl [:] "DecayEn.d" u 1:cc w lp lw 2, "DecayMa.d" u 1:cc w lp lw 2, "DecayFockEn.d" u 1:cc w lp lw 2, "DecayFockMa.d" u 1:cc w lp lw 2, f1(.26,x)
set term x11 3; cc=7; pl [:] "DecayEn.d" u 1:cc w lp lw 2, "DecayMa.d" u 1:cc w lp lw 2, "DecayEn.d" u 1:8 w lp lw 2, "DecayMa.d" u 1:8 w lp lw 2, e2(1.25,x)
set term x11 4; cc= 9; pl [:] "DecayEn.d" u 1:cc w lp lw 2, "DecayMa.d" u 1:cc w lp lw 2, e1(  1,x)
set term x11 5; cc=10; pl [:] "DecayEn.d" u 1:cc w lp lw 2, "DecayMa.d" u 1:cc w lp lw 2, e1(-.5,x)
set term x11 6; cc=7; pl [:] "DecayFockEn.d" u 1:cc w lp lw 2, "DecayFockMa.d" u 1:cc w lp lw 2, "DecayFockEn.d" u 1:8 w lp lw 2, "DecayFockMa.d" u 1:8 w lp lw 2, e2(8,x)


pause -1

