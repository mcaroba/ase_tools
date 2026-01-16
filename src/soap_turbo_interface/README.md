Step 1: Build Fortran library

make

Step 2: Run f90wrap

python3 -m f90wrap src/soap_turbo*.f90 -m soap_turbo -k kind_map.json -v

Step 3: Build Python module with Meson

meson setup builddir --wipe
meson compile -C builddir
meson install -C builddir

Step 4: Try to use soap_turbo in a Python program
python3 test_get_soap.py


Try again:
make clean



Debugging:

Check Python path
python3 -c "import sys; print(sys.path)"

Check Python version:
which python3

Other:
ls -ld /usr/local/lib/python3.13/site-packages
ls -l /usr/local/lib/python3.13/site-packages/soap_turbo
find builddir -name "*.so"
find /usr/local/lib/python3.13/site-packages -name "*.so"
