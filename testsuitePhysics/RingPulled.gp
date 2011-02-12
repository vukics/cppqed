COMP_I={0,1}

complex(x,y)=x*{1,0}+y*COMP_I

conj(x)=complex(real(x),-imag(x))

Unot1=24
Unot2=18

omrec=1e-6

v0=1600

deltaC1=-6
deltaC2=-5

kappa1=60
kappa2=50

eta1={.9,-.2}
eta2={.4,-.4}

position(x)=2*omrec*v0*x

f1(x)=exp( COMP_I*x)
f2(x)=exp(-COMP_I*x)

A=sqrt(Unot1*Unot2)

g(x)=A*conj(f1(x))*f2(x)


Z1(x)=complex(-kappa1,deltaC1-Unot1*abs(f1(x))**2)
Z2(x)=complex(-kappa2,deltaC2-Unot2*abs(f2(x))**2)

det(x)=Z1(x)*Z2(x)+abs(g(x))**2

a1(x)=-(Z2(x)*eta1+COMP_I*     g(x) *eta2)/det(x)
a2(x)=-(Z1(x)*eta2+COMP_I*conj(g(x))*eta1)/det(x)


set term x11 0

plot [.2:][0:] "RingPulled.d" u 1:5, position(x)

set term x11 1

plot [.2:] "RingPulled.d" u 1:9 , "RingPulled.d" u 1:(real(a1($5))) w l

set term x11 2

plot [.2:] "RingPulled.d" u 1:10, "RingPulled.d" u 1:(imag(a1($5))) w l

set term x11 3

plot [.2:] "RingPulled.d" u 1:13, "RingPulled.d" u 1:(real(a2($5))) w l

set term x11 4

plot [.2:] "RingPulled.d" u 1:14, "RingPulled.d" u 1:(imag(a2($5))) w l

pause -1


f1(x)=sin(x)
f2(x)=cos(x)

set term x11 0

plot [.2:][0:] "RingPulledSC.d" u 1:5, position(x)

set term x11 1

plot [.2:] "RingPulledSC.d" u 1:9 , "RingPulledSC.d" u 1:(real(a1($5))) w l

set term x11 2

plot [.2:] "RingPulledSC.d" u 1:10, "RingPulledSC.d" u 1:(imag(a1($5))) w l

set term x11 3

plot [.2:] "RingPulledSC.d" u 1:13, "RingPulledSC.d" u 1:(real(a2($5))) w l

set term x11 4

plot [.2:] "RingPulledSC.d" u 1:14, "RingPulledSC.d" u 1:(imag(a2($5))) w l

pause -1



f1(x)=sin(  x)
f2(x)=sin(2*x)

set term x11 0

plot [.2:][0:] "RingPulledS2S.d" u 1:5, position(x)

set term x11 1

plot [.2:] "RingPulledS2S.d" u 1:9 , "RingPulledS2S.d" u 1:(real(a1($5))) w l

set term x11 2

plot [.2:] "RingPulledS2S.d" u 1:10, "RingPulledS2S.d" u 1:(imag(a1($5))) w l

set term x11 3

plot [.2:] "RingPulledS2S.d" u 1:13, "RingPulledS2S.d" u 1:(real(a2($5))) w l

set term x11 4

plot [.2:] "RingPulledS2S.d" u 1:14, "RingPulledS2S.d" u 1:(imag(a2($5))) w l

pause -1

