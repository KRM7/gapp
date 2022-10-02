#!/bin/sh
  
git clone https://github.com/catchorg/Catch2.git
cd Catch2
git checkout tags/v3.1.0
cmake -Bbuild -H. -DBUILD_TESTING=OFF
sudo env "PATH=$PATH" cmake --build build/ --target install
cd ..