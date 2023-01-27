# ase_tools

Fast Fortran libraries to be used in combination with ASE

**ase_tools** is copyright (c) 2022-2023 by Rina Ibragimova and Miguel A. Caro. It is
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

Add the source directory to your Python path to use the libraries directly:

    echo "export PYTHONPATH=$(pwd)/src:\$PYTHONPATH" >> ~/.bashrc
    source ~/.bashrc

## Available tools

### Identifying surface atoms

Assuming you have an ASE `atoms = Atoms(...)` object:

    from ase_tools import surface_list
    surf_list = surface_list(atoms, r_min, r_max, n_tries, cluster=True)

`surface_list()` returns a list of indices for those atoms identified as surface atoms by the rolling-sphere
algorithm. `cluster=True` is currently required since the implementation does not (yet) support periodic
boundary conditions. `n_tries` is the number of probes, which are placed randomly within the simulation
box. `r_min` and `r_max` are defined graphically below. `r_max` should be slightly bigger than `r_min`.
The higher the value of `n_tries` the more accurate the estimate (and the closer `r_max` can be to
`r_min`).

![Rooling-sphere algorithm](docs/img/rolling_sphere.png)
