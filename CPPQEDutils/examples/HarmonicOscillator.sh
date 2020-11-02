$1/HarmonicOscillatorExact       --T -1 --G 1.          > g1.d
$1/HarmonicOscillatorExact       --T -1 --G  .99        > g1me.d
$1/HarmonicOscillatorExact       --T -1 --G 1.01        > g1pe.d
$1/HarmonicOscillatorRandomTimes --T -1 --G 1.01        > rt.d
$1/HarmonicOscillatorComplex     --T -1 --G 1.01 --dc 1 > co.d
# The gnuplot command---pl [:15] 'g1.d' u 1:3 '%lf %lf (%lf,%lf) (%lf,%lf)' w l, 'g1pe.d' u 1:3 '%lf %lf (%lf,%lf) (%lf,%lf)' w l, 'g1me.d' u 1:3 '%lf %lf (%lf,%lf) (%lf,%lf)' w l, 'rt.d' u 1:3 '%lf %lf (%lf,%lf) (%lf,%lf)', 'co.d' u 1:3 '%lf %lf (%lf,%lf) (%lf,%lf)' w lp