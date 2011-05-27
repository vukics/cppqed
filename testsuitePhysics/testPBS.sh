echo "#!/bin/bash" > temp.sh
echo "#" >> temp.sh
echo "#PBS -N quantologia_testsuitePhysics" >> temp.sh
echo "#PBS -l nodes=1,cput=10:00:00" >> temp.sh
echo "#PBS -q defroute" >> temp.sh

echo "$PATH_TO_EXECS/$2" >> temp.sh
qsub -o $1 temp.sh
#rm temp.sh
