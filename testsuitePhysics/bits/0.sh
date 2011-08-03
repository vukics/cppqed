###############################################
# Two-level atom with different implementations
###############################################

$1 "PTLA_Ev.d" "PTLA_Evolved $ARGS0              "
$1 "PTLA_Ma.d" "PTLA_C++QED  $ARGS0 --evol master"

$1 "PTLA_En.d"   "PTLA_C++QED     $ARGS0               $ARGSEnsemble"

$1 "PLQ_En.d"    "PumpedLossyQbit $ARGS0               $ARGSEnsemble"
$1 "PLQ_EnSch.d" "PumpedLossyQbit $ARGS0 --picture Sch $ARGSEnsemble"
$1 "PLQ_EnUIP.d" "PumpedLossyQbit $ARGS0 --picture UIP $ARGSEnsemble"
$1 "PLQ_MaSch.d" "PumpedLossyQbit $ARGS0 --picture Sch --evol master"
$1 "PLQ_MaUIP.d" "PumpedLossyQbit $ARGS0 --picture UIP --evol master"
