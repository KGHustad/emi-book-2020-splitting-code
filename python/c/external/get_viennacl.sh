#!/bin/bash

set -o errexit

URL=https://downloads.sourceforge.net/project/viennacl/1.7.x/ViennaCL-1.7.1.tar.gz
if [ "$(uname -s)" = "Linux" ]; then
    wget $URL
elif [ "$(uname -s)" = "Darwin" ]; then
    curl -OL $URL
else
    echo "Unsupported platform"
    exit 1
fi

tar xvf ViennaCL-1.7.1.tar.gz
mv ViennaCL-1.7.1/viennacl .
rm -rf ViennaCL-1.7.1
