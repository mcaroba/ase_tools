#!/bin/bash

rm src/*.so

cd src

python3 -m numpy.f2py --f90flags='-fopenmp' -lgomp -c cluster.f90 -m cluster
python3 -m numpy.f2py --f90flags='-fopenmp' -lgomp -c surface.f90 -m surface


# Build soap_turbo interface here conditionally on having cloned the soap_turbo
# submodule (git pull --recurse-submodules)
if [[ -d soap_turbo ]]; then

set -e

PYTHON="${PYTHON:-python3}"
BUILDDIR="builddir"

if [ "$1" = "--clean" ]; then
    echo "Cleaning build artifacts..."
    rm -rf "$BUILDDIR"
    rm -f f90wrap_*.f90 soap_turbo.py libsoap_turbo.a src/*.o src/*.mod
fi

# Detect the prefix of the currently active Python (works with system, venv, Conda)
PREFIX=$("$PYTHON" -c "import sys; print(sys.prefix)")
echo "Installing into Python environment: $PREFIX"
echo "Python executable: $("$PYTHON" --version)"

meson setup "$BUILDDIR" --wipe --prefix "$PREFIX"
meson compile -C "$BUILDDIR"
meson install -C "$BUILDDIR"

# Move the Python module to the parent directory
cd ..
mv src/ase_tools_for_soap.py .

fi

