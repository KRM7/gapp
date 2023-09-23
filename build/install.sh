#!/bin/sh

echo -e "Installing gapp...\n"

cmake --version

BUILD_DIR=$(dirname $(realpath "$0"))
echo -e "\nThe build directory is ${BUILD_DIR}.\n"

echo -e "Installing Debug configuration.\n"
cmake -S $BUILD_DIR/.. -B $BUILD_DIR -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=OFF "$@"
cmake --build $BUILD_DIR --config Debug --parallel 8
cmake --install $BUILD_DIR --config Debug

echo -e "Installing RelWithDebInfo configuration.\n"
cmake -S $BUILD_DIR/.. -B $BUILD_DIR -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=OFF "$@"
cmake --build $BUILD_DIR --config RelWithDebInfo --parallel 8
cmake --install $BUILD_DIR --config RelWithDebInfo

echo -e "Installing Release configuration.\n"
cmake -S $BUILD_DIR/.. -B $BUILD_DIR -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF "$@"
cmake --build $BUILD_DIR --config Release --parallel 8
cmake --install $BUILD_DIR --config Release
