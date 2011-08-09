echo "$2" > temp.sh
qsub -cwd -q all.q -o $1 ./temp.sh
#-N "$1" 
rm temp.sh
