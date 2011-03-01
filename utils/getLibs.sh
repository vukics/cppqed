if [ $1 == "local" ]; then
  echo "Installing locally"
else
  echo "Installing globally"
fi


SCRIPTDIR=`pwd`

cvs -d:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz login
cvs -z3 -d:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz co -P blitz
cd blitz
autoreconf -vif
if [ $1 == "local" ]; then
  echo "Installing Blitz++ locally"
  ./configure --with-pic --prefix=$PWD
  make lib
  make install
else
  echo "Installing Blitz++ globally"
  ./configure --with-pic
  make lib
  sudo make install
fi
cd $SCRIPTDIR

cvs -d:pserver:anonymous@flens.cvs.sourceforge.net:/cvsroot/flens login
cvs -z3 -d:pserver:anonymous@flens.cvs.sourceforge.net:/cvsroot/flens co -P FLENS-lite
cd FLENS-lite
if [ $1 == "local" ]; then
  echo "Installing FLENS locally"
  echo "PREFIX="$PWD >> config
  cat config.ubuntu >> config
else
  echo "Installing FLENS globally"
  cp config.ubuntu config
fi
sed -i -e 's/LDFLAGS    += -llapack -latlas/LDFLAGS    += -llapack/g' -e 's/CXXFLAGS   +=/CXXFLAGS   += -DGSL_CBLAS/g' config
sed -i -e 's_$(PWD)_'$PWD'_g' -e 's_SUBDIRS = $(LIBDIRS) $(EXECDIRS)_SUBDIRS = $(LIBDIRS)_g' Makefile.common
make
if [ $1 == "local" ]; then
  make install
else
  sudo make install
fi
cd $SCRIPTDIR




if [ $1 == "local" ]; then
  cp Jamfile Jamfile.global
  echo "project : usage-requirements <include>blitz/include <include>FLENS-lite/include ;" > Jamfile
  echo "lib blitzlib : : <file>"$SCRIPTDIR"/blitz/lib/libblitz.a ;" >> Jamfile
  echo "lib flenslib : : <file>"$SCRIPTDIR"/FLENS-lite/lib/libflens.so <link>shared ;" >> Jamfile
fi
