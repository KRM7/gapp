#!/bin/bash

echo -e "Installing Catch2...\n"

cmake --version

BUILD_DIR=$(dirname $(realpath "$0"))/../build
echo -e "\nThe build directory is ${BUILD_DIR}.\n"

git clone -b v3.3.0 https://github.com/catchorg/Catch2.git $BUILD_DIR/Catch2
CATCH_BUILD_DIR=$BUILD_DIR/Catch2/build
echo -e "\nThe Catch build directory is ${CATCH_BUILD_DIR}.\n"

echo -e "Installing $1 configuration.\n"
cmake -S $CATCH_BUILD_DIR/.. -B $CATCH_BUILD_DIR -DCMAKE_BUILD_TYPE=$1 -DBUILD_TESTING=OFF "${@:2}"
cmake --build $CATCH_BUILD_DIR --config $1 --parallel
cmake --install $CATCH_BUILD_DIR --config $1
