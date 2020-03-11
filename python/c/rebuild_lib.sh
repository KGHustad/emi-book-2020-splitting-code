#!/usr/bin/env bash

set -o errexit

rm -rf build &&
mkdir -p build &&
pushd build &&
cmake .. &&
make -j4 &&
popd
