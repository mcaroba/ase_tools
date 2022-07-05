# ase_tools

Fast Fortran libraries to be used in combination with ASE

**ase_tools** is copyright (c) 2022 by Rina Ibragimova and Miguel A. Caro. It is
distributed under the GNU General Public License version 3. It relies on a working
installation of ASE and Numpy.

See the LICENSE.md file for detailed information on this
software's license.

## Installation

### Prerrequisites

- Numpy
- A Fortran compiler (successfully tested with `gfortran`)

### Building the libraries

Clone the **ase_tools** repository:

    git clone http://github.com/mcaroba/ase_tools.git

Execute the build script:

    cd ase_tools/
    ./build_libraries.sh

Add the root directory to your Python path:

    echo "export PYTHONPATH=$(pwd):\$PYTHONPATH" >> ~/.bashrc
    source ~/.bashrc
