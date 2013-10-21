for i in `seq 1 1 3`; do time bjam -ad0 C++QED toolset=gcc-4.4 release; done
for i in `seq 1 1 3`; do time bjam -ad0 1particle1mode.o toolset=gcc-4.4 release; done
for i in `seq 1 1 3`; do time bjam -ad0 4qbits.o toolset=gcc-4.4 release; done
bjam -j7 4qbits 1particle1mode toolset=gcc-4.4 release
for i in `seq 1 1 3`; do time bin/scripts/gcc-4.4/release/1particle1mode >/dev/null; done
for i in `seq 1 1 3`; do time bin/scripts/gcc-4.4/release/4qbits --eta 1 >/dev/null; done
