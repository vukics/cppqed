###########################
# Non-interacting composite
###########################

$1 "QMJ_En.d"    "QbitMode_C++QED  $ARGS3 $ARGSEnsemble"
$1 "QMJ_Ma.d"    "QbitMode_C++QED  $ARGS3 --evol master"
# This is an example where UIP is slower than Sch because the largest frequency is the pump

$1 "QMJ_MaSch.d" "QbitMode_C++QED  $ARGS3 --evol master --picture Sch"

