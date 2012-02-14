for evalScript in Decay.py PTLAHL.py PTLA.py PumpedLossyMode.py QMJ_Int.py PLM_FT.py freeParticle.py ; do pyEvalScripts/$evalScript $1 $2 2>/dev/null ;
# redirect stderr to /dev/null
#| awk '{if (index($0,"maxiter")==0) print}'
done
