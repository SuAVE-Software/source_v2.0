#!/bin/bash


# Installing LIBXDR Library ======================
cd library/libxdrfile-master/
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
sudo make install 
cd ../../..
# ================================================


# Installing GMXXDR Library ======================
cd library/libgmxfort-master/
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
sudo make install
cd ../../..
# ================================================



