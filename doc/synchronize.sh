figures
bash processImages.sh
cd ..

make html
make latexpdf

rsync -Cavuz --exclude '*~' --delete _build/html/* vukics,cppqed@web.sourceforge.net:htdocs/
