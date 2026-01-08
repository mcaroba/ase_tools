# ase_tools is copyright (c) 2022-2026 of Miguel Caro and Rina Ibragimova
# See License.md for more information

import numpy as np
from ase import Atoms
from src.surface import surface_tools
from src.cluster import cluster_module
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write,read
import subprocess

###############################################################################################
#
# This function takes an ASE atoms object as argument and returns a list of surface atoms
# according to a probe-sphere algorithm based on Monte Carlo sampling (test spheres are
# randomly placed in the vacuum side of the simulation box).
#
# r_min  :: this is the minimum distance allowed between an atom and a probe sphere in Ang.
# r_max  :: this is the maximum distance allowed between an atom and a probe sphere in Ang.
# n_tries :: this is the number of probe spheres that will be used to probe the simulation
#           box
#
def surface_list(atoms, r_min, r_max, n_tries, cluster=False):
    if not cluster:
        raise ValueError("Currently the surface_list function can only be used for cluster systems!")
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    is_surface = surface_tools.surface_atoms(pos, cell, r_min, r_max, n_tries, cluster)
    surf_list = []
    for i in range(0, len(atoms)):
        if is_surface[i] == 1:
            surf_list.append(i)
    return surf_list
###############################################################################################
#
# This function splits a series of atoms into individual molecules, according to some bonding
# cutoffs, and returns a list of Atoms() objects each containin a single molecule.
#
def split_atoms(atoms, bonding_cutoff):
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    symb = atoms.get_chemical_symbols()
#   Assert this is an orhthorhombic unit cell - THIS SHOULD BE EXTENDED TO TRICLINIC!!!
    lx = cell[0][0]
    ly = cell[1][1]
    lz = cell[2][2]
    if not (cell[0][1] == 0 and cell[0][2] == 0 and cell[1][0] == 0 and \
            cell[1][2] == 0 and cell[0][1] == 0 and cell[0][2] == 0):
        raise ValueError("Currently split_atoms only works with orthogonal simulation boxes!")
#   Check whether the bonding cutoff is a scalar, an array, or a dictionary
    if isinstance(bonding_cutoff, float) or isinstance(bonding_cutoff, int):
        cutoffs = np.zeros(len(atoms)) + bonding_cutoff
    elif len(bonding_cutoff) == len(atoms) and all([(isinstance(a, float) or isinstance(a, int)) for a in bonding_cutoff]):
        cutoffs = bonding_cutoff
    elif type(bonding_cutoff) is dict:
        cutoffs = np.zeros(len(atoms))
        for i in range(0, len(atoms)):
            this_symb = symb[i]
            try:
                cutoffs[i] = bonding_cutoff[this_symb]
            except:
                raise ValueError("Cutoff not defined for species " + this_symb)
    else:
        raise ValueError("bonding_cutoff must be a scalar, an array with the same dimensions as the number of atoms, or a dictionary!")
#   Build the list of indices
    list = cluster_module.cluster_atoms(np.transpose(pos), lx, ly, lz, cutoffs)
#   Build the lists of Atoms() objects
    atoms_list = []
    for j in range(1, max(list)+1):
        pos_new = []
        symb_new = ""
        for i in range(0, len(list)):
            if list[i] == j:
                pos_new.append(pos[i])
                symb_new += symb[i]
        atoms_list.append(Atoms(symb_new, cell = cell, positions=pos_new, pbc=True))
#   Return the list of Atoms() objects
    return atoms_list
###############################################################################################
#
# This function reads a PDB file and optionally removes spurious periodic copies of the atoms
#
def read_pdb(filename, remove_spurious_atoms=False):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    symbols = []
    positions = []
    for line in lines:
        words = line.split()
        if len(words) > 0 and words[0] == "CRYST1":
            a, b, c, alpha, beta, gamma = [float(num) for num in words[1:7]]
        if len(words) > 0 and words[0] == "ATOM":
            symb = words[2]
            pos = [float(num) for num in words[3:6]]
            symbols.append(symb); positions.append(pos)
    atoms = Atoms(symbols, positions=positions, cell = [a, b, c, alpha, beta, gamma], pbc=True)
    if remove_spurious_atoms:
        for i in range(len(atoms)-1, -1, -1):
            for j in range(i-1, -1, -1):
                d = atoms.get_distance(i, j, mic=True)
                if d < 0.01:
                    del atoms[i]
                    break
    return atoms
###############################################################################################
#
# This function creates an incredibly inefficient interface to the TurboGAP code. Note that it
# requires for the user to have a local input file and local gap_files/ directory. It will return
# the energy and forces of the last image in the trajectory_out.xyz. It only makes sense to use
# it with a minimal input file which you would usually use to run "turbogap predict"
#
class TurboGAP(Calculator):

    implemented_properties = ['energy', 'forces']

    def __init__(self, atoms=None, verbose=False):
        Calculator.__init__(self)
        self.verbose = verbose

    def calculate(self,
                  atoms=None,
                  properties=['energy', 'forces'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        write("atoms.xyz", self.atoms)
        if not self.verbose:
            subprocess.run(["turbogap", "predict"], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        else:
            subprocess.run(["turbogap", "predict"])
        atoms_temp = read("trajectory_out.xyz", index=-1)
        self.results['energy'] = atoms_temp.get_potential_energy()
        self.results['forces'] = atoms_temp.get_forces()
###############################################################################################
