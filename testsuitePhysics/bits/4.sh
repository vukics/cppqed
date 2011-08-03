#################################
# Pure decay in composite systems
#################################

$1 "DecayMa.d"     "QbitMode_C++QED $ARGS4 --minit '(1.,-.5)' --evol master  "
$1 "DecayEn.d"     "QbitMode_C++QED $ARGS4 --minit '(1.,-.5)' $ARGSEnsemble  "

$1 "DecayFockMa.d" "QbitMode_C++QED $ARGS4 --minitFock 8 --evol master  "
$1 "DecayFockEn.d" "QbitMode_C++QED $ARGS4 --minitFock 8 $ARGSEnsemble  "
