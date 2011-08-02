if [ $1 == "local" ]; then
  echo "Installing locally"
else
  echo "Installing globally"
fi


SCRIPTDIR=`pwd`

cvs -d:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz login
cvs -z3 -d:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz co -P blitz
# hg clone http://blitz.hg.sourceforge.net:8000/hgroot/blitz/blitz
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
sed -i '
s|LDFLAGS    += -llapack -latlas|LDFLAGS    += -llapack|g
s|CXXFLAGS   +=|CXXFLAGS   += -DGSL_CBLAS|g' config
sed -i '
s|$(PWD)|'$PWD'|g
s|SUBDIRS = $(LIBDIRS) $(EXECDIRS)|SUBDIRS = $(LIBDIRS)|g' Makefile.common
make
if [ $1 == "local" ]; then
  make install
else
  sudo make install
fi
cd $SCRIPTDIR




if [ $1 == "local" ]; then
  cp Jamroot Jamroot.global
  sed -i '
s|usage-requirements|usage-requirements <include>blitz/include <include>FLENS-lite/include|g
s|lib flenslib : : <name>flens|lib flenslib : : <file>'$SCRIPTDIR'/FLENS-lite/lib/libflens.so|g
s|lib blitzlib : : <name>blitz|lib blitzlib : : <file>'$SCRIPTDIR'/blitz/lib/libblitz.a|g' Jamroot
fi
