# exit on first error
set -e

if [ "$1" == "local" ]; then
  echo "Installing locally"
else
  echo "Installing globally"
fi


SCRIPTDIR=`pwd`
LOCALLIBDIR=$SCRIPTDIR/local

hg clone http://blitz.hg.sourceforge.net:8000/hgroot/blitz/blitz
cd blitz
autoreconf -vif
if [ "$1" == "local" ]; then
  echo "Installing Blitz++ locally"
  ./configure --with-pic --enable-serialization --prefix=$LOCALLIBDIR
  make lib
  make install
else
  echo "Installing Blitz++ globally"
  ./configure --with-pic --enable-serialization
  make lib
  sudo make install
fi
cd $SCRIPTDIR

cvs -z3 -d:pserver:anonymous:@flens.cvs.sourceforge.net:/cvsroot/flens co -P FLENS-lite
cd FLENS-lite
if [ "$1" == "local" ]; then
  echo "Installing FLENS locally"
  echo "PREFIX="$LOCALLIBDIR >> config
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
if [ "$1" == "local" ]; then
  make install
else
  sudo make install
fi
cd $SCRIPTDIR


echo -e "\n\nInstallation of blitz++ and FLENS successfull.\n\n"

if [ "$1" == "local" ]; then
  echo -e "****************************************************************************************************"
  echo -e "Use the following commands to configure the build system (build type can be 'debug' or 'release'):\n"
  echo -e "mkdir -p build; cd build; cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_PREFIX_PATH=$LOCALLIBDIR .."
  echo -e "****************************************************************************************************\n"

  echo -e "****************************************************************************************************"
  echo -e "For Boost.Build, add the following lines to ~/user-config.jam\n"
  echo -e "local PREFIX = $LOCALLIBDIR ;"
  echo -e 'project genconfig : requirements <include>$(PREFIX)/include <library-path>$(PREFIX)/lib ;'
  echo -e "****************************************************************************************************\n"
fi
