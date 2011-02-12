for dir in tutorial manual; do
  cd $dir/figures
  bash ../../processImages.sh
  cd ../..
done

cd tutorial; make latexpdf; 
cd ..

make html

rsync -Cavuz --exclude '*~' --delete _build/html/* vukics,cppqed@web.sourceforge.net:htdocs/
