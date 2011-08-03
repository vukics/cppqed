#########################################
# Two-level atom with Heisenberg-Langevin
#########################################

$1 "PTLAHL_Si.d" "PTLA_C++QED    $ARGS1 --no_noise"
$1 "PTLAHL_Ev.d" "PTLA_EvolvedHL $ARGS1           "
