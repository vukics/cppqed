###########################################
# Lossy mode with different implementations
###########################################

$1 "PLM_Ev.d"    "PumpedLossyMode_Evolved $ARGS2                            "
$1 "PLM_En.d"    "PumpedLossyMode_C++QED  $ARGS2               $ARGSEnsemble"
$1 "PLM_Si.d"    "PumpedLossyMode_C++QED  $ARGS2                            "
$1 "PLM_SiSch.d" "PumpedLossyMode_C++QED  $ARGS2 --picture Sch              "
$1 "PLM_SiUIP.d" "PumpedLossyMode_C++QED  $ARGS2 --picture UIP              "
$1 "PLM_MaSch.d" "PumpedLossyMode_C++QED  $ARGS2 --picture Sch --evol master"
$1 "PLM_MaUIP.d" "PumpedLossyMode_C++QED  $ARGS2 --picture UIP --evol master"
