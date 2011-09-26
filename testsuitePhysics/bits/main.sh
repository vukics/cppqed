export pathToExecs=$2

export pluginToRun=$3

export dataDirectory=$4
export additionalArgs=$5 

# Shared variables:

export ARGSTraj="--dc 0 --Dt .01 --T 2"
# The fifth positional parameter being an additional set of parameters

export ARGSEnsemble="--nTraj 200 --evol ensemble"

export ARGS0Base="--etat '(1,3)' --gamma 1.2 --deltaA -1 --qbitInit '(.4,.3)'"

export ARGS2Base="--minit '(1,-.5)' --eta '(30,10)' --deltaC 100 --cutoff 30"

# The l function:

l ()
{ 
bash $pluginToRun "$dataDirectory/$1" "$pathToExecs/$2 $additionalArgs"
}

export -f l

# Calling the main script:

bash $1

