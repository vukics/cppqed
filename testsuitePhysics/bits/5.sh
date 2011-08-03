####################
# Interacting binary
####################


$1 "QMJ_Int_Ev.d"     "QbitMode_Evolved $ARGS3                         "

$1 "QMJ_Int_Matrix.d" "QbitMode_Matrix  $ARGS3               --no_noise"

$1 "QMJ_Int_Si.d"     "QbitMode_C++QED  $ARGS3               --no_noise"

$1 "QMJ_Int_En.d"     "QbitMode_C++QED  $ARGS3 $ARGSEnsemble"
$1 "QMJ_Int_Ma.d"     "QbitMode_C++QED  $ARGS3 --evol master"

$1 "QMJ_Int_MaSch.d"  "QbitMode_C++QED  $ARGS3 --evol master --picture Sch"

$1 "QMJ_Int_SiSch.d"  "QbitMode_C++QED  $ARGS3 --picture Sch --no_noise"
$1 "QMJ_Int_SiUIP.d"  "QbitMode_C++QED  $ARGS3 --picture UIP --no_noise"
