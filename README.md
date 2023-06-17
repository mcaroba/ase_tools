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

    echo "export PYTHONPATH=$(pwd):\$PYTHONPATH" >> ~/.bashrc
    source ~/.bashrc

## ASE tools

These are the general tools compatible with a general ASE `Atoms()` object.

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

NOTE: currently, this tool does not allow you to handle systems under periodic boundary conditions, i.e.,
you'll need to pass the `cluster=True` option as argument, otherwise the code will throw an error.

### Splitting a system into individual molecules

If you have a simulation box with a series of atoms bonded such that they form molecules, and want to
have those molecules separated into individual ASE's `Atoms(...)` objects, you can use ase_tools'
`split_atoms(...)` function:

    from ase_tools import split_atoms
    db = split_atoms(atoms, bonding_cutoff={"H": 1.3, "C": 1.8})

`split_atoms()` returns a list of `Atoms()` objects, each with the molecules that could be constructed
by assuming two atoms are bonded if `distance[i,j] < (cutoff[i]+cutoff[j])/2.`. The bonding cutoff can
be a scalar, an array with the same length as the number of atoms in the input `Atoms()` object, or
a dictionary containing a cutoff value for each species present in the system, as in the example above.

## ASE tools for VASP

These are tools which are used to facilitate calculations with VASP.

### Generating KPOINTS files accounting for the presence of vacuum

This function makes a VASP `KPOINTS` file compatible with the chosen combination of `KSPACING`
and `KGAMMA` parameters. You can also apply a shift to the k-mesh origin. The main
functionality is that it can detect the presence of vacuum in the simulation cell, such
that the sampling along the vacuum direction(s) will be one k point only. This can help
to save CPU time in situations where a high-throughput approach is being used and the user
wants to avoid unnecessarily dense k-point sampling for surfaces, nanoparticles, etc.

The `vacuum_check` variable enables detection of vacuum, and the `vacuum_cutoff` variable tells
the code what gap in the atomic density corresponds to the presence of vacuum. E.g.,
`vacuum_cutoff = 5`. means that whenever there is a 5 Angstrom gap with no atoms present
along some spatial direction, the corresponding k-sampling will be set to 1 along that
direction. When this leads to Gamma-point sampling (e.g., for nanoparticles and molecules),
the user should further use this information to use the Gamma-point VASP binary, for extra
CPU time savings.

Currently, the code can only detect trivial vacuum gaps perpendicular to a given plane (where
this plane is defined by two of the lattice vectors through their cross product). Complex
"curved" gaps cannot be detected. This should not be an issue for most applications.

Example usage, where `atoms` is a valid ASE `Atoms()` object:

    from ase_tools_for_vasp import make_kpoints_file
    make_kpoints_file(atoms, filename="KPOINTS", kspacing=0.5, kgamma=True,
                      shift=[0.,0.,0.], vacuum_check=True, vacuum_cutoff=5.):
