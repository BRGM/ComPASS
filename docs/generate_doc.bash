#!/bin/bash

rootdir=$1
wheeltag=$2

/bin/bash sdk/install_wheel.bash ${wheeltag}

pushd ${rootdir}/docs/doxygen
doxygen Doxyfile
popd


pushd ${rootdir}/docs
for f in README INSTALL LICENSE
do
    cp -vf ${rootdir}/$f.rst .
done
sphinx-apidoc ../ComPASS -o python_reference
sphinx-build . html
popd
