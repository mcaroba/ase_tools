#!/bin/bash

rm src/*.so

cd src

python3 -m numpy.f2py --f90flags='-fopenmp' -lgomp -c cluster.f90 -m cluster
python3 -m numpy.f2py --f90flags='-fopenmp' -lgomp -c surface.f90 -m surface


# Build soap_turbo interface here conditionally on having cloned the soap_turbo
# submodule
if [[ -d src/soap_turbo ]]; then

fi
