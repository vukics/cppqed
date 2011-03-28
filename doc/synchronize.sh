cd figures
bash processImages.sh
cd ..

make latexpdf
make html

rsync -Cavuz --exclude '*~' --delete _build/html/* vukics,cppqed@web.sourceforge.net:htdocs/
