echo $1
echo "time $2 > $1" > temp.sh
chmod u+x temp.sh
./temp.sh
rm temp.sh
