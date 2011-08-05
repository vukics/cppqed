echo "#!/bin/bash" > temp.sh
echo "#" >> temp.sh
echo "#PBS -N $1" >> temp.sh
#echo "#PBS -l nodes=1,cput=10:00:00" >> temp.sh
echo "#PBS -q eternity -e errorfiles" >> temp.sh

echo "$2 --o $1" >> temp.sh
qsub -o /dev/null temp.sh
rm temp.sh
