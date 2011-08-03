echo "$PATH_TO_EXECS/$2" > temp.sh
qsub -cwd -q all.q -o $1 -N $1 ./temp.sh
rm temp.sh
