# Starting from fresh Ubuntu 20.04 LTS installation
# Install prerequisites:
sudo apt install libgsl-dev libboost-all-dev libeigen3-dev clang cmake cmake-curses-gui git python3-scipy # ipython3
# Clone repo + submodules:
git clone --recurse-submodules --remote-submodules https://github.com/vukics/cppqed.git
cd cppqed
# Checkout Ubuntu branch: 
git checkout Ubuntu20.04LTS
mkdir build; cd build
# Configure build:
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
# Build everything needed for tests:
make -j4 boostTester fewer_scripts compileTests cpypyqed 
# Tests can be run with
# ctest -j4
# from the main build folder
