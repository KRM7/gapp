#!/bin/sh

git clone https://github.com/fmtlib/fmt.git
cd fmt
git checkout tags/9.0.0
cmake -Bbuild -S. -DFMT_TEST=OFF
cmake --build build/ --target install --config Debug
cmake --build build/ --target install --config Release
cmake --build build/ --target install --config RelWithDebInfo
cd ..