echo $1 "          " $2
echo "$2 > $1" > temp.sh
chmod u+x temp.sh
./temp.sh
rm temp.sh
