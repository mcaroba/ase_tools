#!/bin/bash

rm src/*.so

cd src

python3 -m numpy.f2py --f90flags='-fopenmp' -lgomp -c cluster.f90 -m cluster


