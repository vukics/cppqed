export PATH_TO_EXECS=/home/vukics/quantologia/bin/scripts/gcc-4.5/release
#export PATH_TO_EXECS=/nfshome/vukics/h/quantologia/bin/scripts/gcc-4.4/release

PATH=.:$PATH_TO_EXECS:"${PATH}"

ARGSTrajBase=""
ARGSTraj="$ARGSTrajBase --dc 0 --Dt .01 --T 2"

ARGSEnsemble="--nTraj 200 --evol ensemble"

ARGS0Base="--etat '(1,3)' --gamma 1.2 --deltaA -1 --qbitInit '(.4,.3)'"
ARGS0="$ARGSTraj $ARGS0Base"

ARGS1Base="--etat '(5,-4)' --gamma 100 --deltaA 1000 --qbitInit '(-.1,.05)'"
ARGS1="$ARGSTraj $ARGS1Base --Dt .001"

ARGS2Base="--minit '(1,-.5)' --eta '(30,10)' --deltaC 100 --cutoff 30"
ARGS2="$ARGSTraj $ARGS2Base --T 1"

ARGS3="$ARGSTraj $ARGS0Base $ARGS2Base --g 0"

ARGS4="$ARGSTraj --qbitInit '(.4,.3)' --etat 0 --gamma 1 --deltaA 0 --eta 0 --kappa 2 --deltaC 0 --g 0"

ARGS5Base="--deltaA -1000 --gamma 100 --g '(80,10)' --etat '(-8,1)' --eta '(1,8)' --kappa 10 --deltaC -62"
ARGS5="$ARGSTraj $ARGS5Base --Dt .001 --T 1"
