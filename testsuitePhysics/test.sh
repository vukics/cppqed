# export PATH_TO_EXECS=/home/vukics/quantologia/bin/scripts/gcc-4.4/release
export PATH_TO_EXECS=/nfshome/vukics/h/quantologia/bin/scripts/gcc-4.4/release

PATH=.:$PATH_TO_EXECS:"${PATH}"

ARGSTraj="--dc 0 --Dt .01 --T 2"

ARGSEnsemble="--nTraj 200 --evol ensemble"

##################################################
# I. Two-level atom with different implementations
##################################################

ARGS0Base="--etat '(1,3)' --gamma 1.2 --deltaA -1 --qbitInit '(.4,.3)'"

ARGS0="$ARGSTraj $ARGS0Base"

$1 "PTLA_Ev.d" "PTLA_Evolved $ARGS0              "
$1 "PTLA_Ma.d" "PTLA_C++QED  $ARGS0 --evol master"

$1 "PTLA_En.d"   "PTLA_C++QED     $ARGS0               $ARGSEnsemble"

$1 "PLQ_En.d"    "PumpedLossyQbit $ARGS0               $ARGSEnsemble"
$1 "PLQ_EnSch.d" "PumpedLossyQbit $ARGS0 --picture Sch $ARGSEnsemble"
$1 "PLQ_EnUIP.d" "PumpedLossyQbit $ARGS0 --picture UIP $ARGSEnsemble"
$1 "PLQ_MaSch.d" "PumpedLossyQbit $ARGS0 --picture Sch --evol master"
$1 "PLQ_MaUIP.d" "PumpedLossyQbit $ARGS0 --picture UIP --evol master"


#############################################
# II. Two-level atom with Heisenberg-Langevin
#############################################

ARGS0bisBase="--etat '(5,-4)' --gamma 100 --deltaA 1000 --qbitInit '(-.1,.05)'"

ARGS0bis="$ARGSTraj $ARGS0bisBase --Dt .001"

$1 "PTLAHL_Si.d" "PTLA_C++QED    $ARGS0bis --no_noise"
$1 "PTLAHL_Ev.d" "PTLA_EvolvedHL $ARGS0bis           "


################################################
# III. Lossy mode with different implementations
################################################

ARGS1Base="--minit '(1,-.5)' --eta '(30,10)' --deltaC 100 --cutoff 30"

ARGS1="$ARGSTraj $ARGS1Base --T 1"

$1 "PLM_Ev.d"    "PumpedLossyMode_Evolved $ARGS1                            "
$1 "PLM_En.d"    "PumpedLossyMode_C++QED  $ARGS1               $ARGSEnsemble"
$1 "PLM_Si.d"    "PumpedLossyMode_C++QED  $ARGS1                            "
$1 "PLM_SiSch.d" "PumpedLossyMode_C++QED  $ARGS1 --picture Sch              "
$1 "PLM_SiUIP.d" "PumpedLossyMode_C++QED  $ARGS1 --picture UIP              "
$1 "PLM_MaSch.d" "PumpedLossyMode_C++QED  $ARGS1 --picture Sch --evol master"
$1 "PLM_MaUIP.d" "PumpedLossyMode_C++QED  $ARGS1 --picture UIP --evol master"

###############################
# IV. Non-interacting composite
###############################

ARGS2="$ARGSTraj $ARGS0Base $ARGS1Base --g 0"

$1 "QMJ_En.d"    "QbitMode_C++QED  $ARGS2 $ARGSEnsemble"
$1 "QMJ_Ma.d"    "QbitMode_C++QED  $ARGS2 --evol master"
# This is an example where UIP is slower than Sch because the largest frequency is the pump

$1 "QMJ_MaSch.d" "QbitMode_C++QED  $ARGS2 --evol master --picture Sch"


####################################
# V. Pure decay in composite systems
####################################

ARGSDecay="--qbitInit '(.4,.3)' --etat 0 --gamma 1 --deltaA 0 --eta 0 --kappa 2 --deltaC 0 --g 0 --dc 0 --Dt .01 --T 2"

$1 "DecayMa.d"     "QbitMode_C++QED $ARGSDecay --minit '(1.,-.5)' --evol master  "
$1 "DecayEn.d"     "QbitMode_C++QED $ARGSDecay --minit '(1.,-.5)' $ARGSEnsemble  "

$1 "DecayFockMa.d" "QbitMode_C++QED $ARGSDecay --minitFock 8 --evol master  "
$1 "DecayFockEn.d" "QbitMode_C++QED $ARGSDecay --minitFock 8 $ARGSEnsemble  "

########################
# VI. Interacting binary
########################

ARGS3Base="--deltaA -1000 --gamma 100 --g '(80,10)' --etat '(-8,1)' --eta '(1,8)' --kappa 10 --deltaC -62"

ARGS3="$ARGSTraj $ARGS3Base --Dt .001 --T 1"

$1 "QMJ_Int_Ev.d"     "QbitMode_Evolved $ARGS3                         "

$1 "QMJ_Int_Matrix.d" "QbitMode_Matrix  $ARGS3               --no_noise"

$1 "QMJ_Int_Si.d"     "QbitMode_C++QED  $ARGS3               --no_noise"

$1 "QMJ_Int_En.d"     "QbitMode_C++QED  $ARGS3 $ARGSEnsemble"
$1 "QMJ_Int_Ma.d"     "QbitMode_C++QED  $ARGS3 --evol master"

$1 "QMJ_Int_MaSch.d"  "QbitMode_C++QED  $ARGS3 --evol master --picture Sch"

$1 "QMJ_Int_SiSch.d"  "QbitMode_C++QED  $ARGS3 --picture Sch --no_noise"
$1 "QMJ_Int_SiUIP.d"  "QbitMode_C++QED  $ARGS3 --picture UIP --no_noise"

####################
# VII. Free particle
####################

ARGS4="--vClass 0 --fin 6 --pinit '(-0.5 10 0.1 0)' --T .5 --dc 0 --Dt 0.01"

$1 "freeParticle.d"     "SingleParticle $ARGS4              "
$1 "freeParticle_Sch.d" "SingleParticle $ARGS4 --picture Sch"
$1 "freeParticle1p1m.d" "1particle1mode $ARGS4 $ARGS1Base --Unot 0 --1p1mconf 4"

###################################################
# VIII. Particle oscillating in classical potential
###################################################

ARGS5="--pinit '(-0.01 0 0 0)' --modePart Cos --fin 10"

$1 "oscillation10000.d" "SingleParticle $ARGS5 --vClass -10000 --T 1 "
$1 "oscillation1000.d"  "SingleParticle $ARGS5 --vClass -1000  --T 3 "
$1 "oscillation100.d"   "SingleParticle $ARGS5 --vClass -100   --T 10"

$1 "oscillation10000_Sch.d" "SingleParticle $ARGS5 --vClass -10000 --T 1 --picture Sch"

#$PATH/SingleParticle $ARGS5 --picture Sch --vClass -100 --T 10 > oscillation100_Sch.d

#$PATH/1particle1mode $ARGS --evol master --fin 8  --etaeff -10000 --T 1 --dc 0 --Dt 0.01 > oscillation10000.dMaster


#####################################################
# IX. Heavy particle not affected by the cavity field
#####################################################

ARGS6Base="--cutoff 3 --fin 11 --deltaC -1e5 --kappa 1e5 --T .1 --pinit '(0 1000 0.05 0)'"

ARGS6="$ARGS6Base --Unot -2e7 --vClass -0.1"

$1 "confAZ_Sin_K1.d" "1particle1mode $ARGS6                 "
$1 "confAZ_Sin_K2.d" "1particle1mode $ARGS6 --kPart 2       "
$1 "confAZ_Cos.d"    "1particle1mode $ARGS6 --modePart Cos  "
$1 "confAZ_Plus.d"   "1particle1mode $ARGS6 --modePart Plus "
$1 "confAZ_Minus.d"  "1particle1mode $ARGS6 --modePart Minus"

ARGS6bis="$ARGS6Base --Unot -6e4 --vClass -10 --1p1mconf 2"

$1 "confAX_Sin_K1.d" "1particle1mode $ARGS6bis                "
$1 "confAX_Sin_K2.d" "1particle1mode $ARGS6bis --kCav 2       "
$1 "confAX_Cos.d"    "1particle1mode $ARGS6bis --modeCav Cos  "
$1 "confAX_Plus.d"   "1particle1mode $ARGS6bis --modeCav Plus "
$1 "confAX_Minus.d"  "1particle1mode $ARGS6bis --modeCav Minus"


######################
# X. Multilevel system
######################

ARGS7="--T 120 --dc 0 --Dt .5"

delta02="-1e3"
reta0="6"
ieta0="5"
reta1="-3"
ieta1="2"
deltaC="-10"
kappa="100"

ARGS7bis="$ARGS7 --gammas '[.01 .005]' --deltas '[0 0 $delta02]'"

$1 "RamanTwoLevel.d"   "PumpedLossyQbit $ARGS7 --gamma 0 --qbitInit '(.7071,0)' `../Raman.py $delta02 $reta0 $ieta0 $reta1 $ieta1`"
$1 "RamanRegression.d" "Raman       $ARGS7bis --etas '[($reta0,$ieta0) ($reta1,$ieta1)]'"
$1 "RamanCavity.d" "CavityRaman $ARGS7bis --etas '[($reta0,$ieta0)]' --deltaC $deltaC --kappa $kappa --cutoff 1200 --Dt .05 --eps 1e-3 --minit '(5,-3)' `../RamanCavity.py $deltaC $kappa $reta1 $ieta1`"

##########
# XI. Ring
##########

ARGS8Base="--pinit '(-.3 0 .01 0)' --UnotM .3 --UnotP .6 --deltaCM -5 --deltaCP -6 --kappaM .5 --kappaP .4 --etaP '(.3,-.1)' --etaM '(.2,-.2)' --T 15 --dc 0 --Dt .1 --fin 10 --cutoffM 3 --cutoffP 3 --omrec 1e-10"

ARGS8="$ARGS8Base --modeCavP Sin --modeCavM Cos"

ARGS8bis="$ARGS8Base --modeCavP Sin --modeCavM Sin --kCavM 2"

$1 "Ring_Ev.d"    "Ring_Evolved $ARGS8Base"
$1 "Ring.d"       "Ring         $ARGS8Base"
$1 "RingSC_Ev.d"  "Ring_Evolved $ARGS8"
$1 "RingSC.d"     "Ring         $ARGS8"
$1 "RingS2S_Ev.d" "Ring_Evolved $ARGS8bis"
$1 "RingS2S.d"    "Ring         $ARGS8bis"


ARGS9Base="--pinit '(0 1600 .05 0)' --UnotM 9 --UnotP 12 --deltaCM -5 --deltaCP -6 --kappaM 50 --kappaP 60 --etaP '(.9,-.2)' --etaM '(.4,-.4)' --T 120 --dc 0 --Dt .1 --fin 12 --cutoffM 3 --cutoffP 3 --omrec 1e-6"

ARGS9="$ARGS9Base --modeCavP Sin --modeCavM Cos"

ARGS9bis="$ARGS9Base --modeCavP Sin --modeCavM Sin --kCavM 2"

$1 "RingPulled.d"    "Ring $ARGS9Base"
$1 "RingPulledSC.d"  "Ring $ARGS9"
$1 "RingPulledS2S.d" "Ring $ARGS9bis"


################
# XII. Tunneling
################

ARGS="--vClass -10 --modePart Sin --T 200 --pinit '(-.5 0 0 0)' --dc 0 --Dt .1"

$1 "tunnel.d"     "SingleParticle $ARGS --picture IP  --fin 8"
$1 "tunnel_Sch.d" "SingleParticle $ARGS --picture Sch --fin 8"


#$PATH/2particles1mode $ARGS --fin 8 > tunnel2particles.d


#########################################
# XIII. NX_CoupledModesElim as oscillator
#########################################

$1 "NXCME_oscillation.d" "NX_coupledModesElim --fin 6 --T 10 --pinit '(-1 0.02 .707 0)' --eta 0"

############################################################################
# XIV. Pure decay in even more composite systems involving alternative decay
############################################################################

ARGSDecay="--eta0 0 --eta1 0 --eta2 0 --eta3 0 --deltaC0 0 --deltaC1 0 --deltaC2 0 --deltaC3 0 --cutoff0 6 --cutoff1 5 --cutoff2 4 --cutoff3 3 --kappa0 1 --kappa1 2 --kappa2 .5 --kappa3 1 --minitFock0 5 --minitFock1 4 --minitFock2 3 --minitFock3 2 --dc 0 --Dt .01 --T 2"

$1 "DecayCompositeMa.d" "FourModes $ARGSDecay --evol master"
$1 "DecayCompositeEn.d" "FourModes $ARGSDecay $ARGSEnsemble"
