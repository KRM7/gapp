#!/bin/bash

echo -e "Installing Catch2...\n"

git clone -b v3.3.0 https://github.com/catchorg/Catch2.git
mkdir Catch2/build

cmake --version

BUILD_DIR=$(dirname $(realpath "$0"))/Catch2/build
echo -e "\nThe build directory is ${BUILD_DIR}.\n"

echo -e "Installing Debug configuration.\n"
cmake -S $BUILD_DIR/.. -B $BUILD_DIR -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=OFF "$@"
cmake --build $BUILD_DIR --config Debug --parallel
cmake --install $BUILD_DIR --config Debug

echo -e "Installing RelWithDebInfo configuration.\n"
cmake -S $BUILD_DIR/.. -B $BUILD_DIR -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=OFF "$@"
cmake --build $BUILD_DIR --config RelWithDebInfo --parallel
cmake --install $BUILD_DIR --config RelWithDebInfo

echo -e "Installing Release configuration.\n"
cmake -S $BUILD_DIR/.. -B $BUILD_DIR -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF "$@"
cmake --build $BUILD_DIR --config Release --parallel
cmake --install $BUILD_DIR --config Release
