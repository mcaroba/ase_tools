# ase_tools is copyright (c) 2022 of Miguel Caro and Rina Ibragimova
# See License.md for more information

import numpy as np
from src.surface import surface_tools

###############################################################################################
#
# This function takes an ASE atoms object as argument and returns a list of surface atoms
# according to a probe-sphere algorithm based on Monte Carlo sampling (test spheres are
# randomly placed in the vacuum side of the simulation box).
#
# r_min  :: this is the minimum distance allowed between an atom and a probe sphere in Ang.
# r_max  :: this is the maximum distance allowed between an atom and a probe sphere in Ang.
# n_tris :: this is the number of probe spheres that will be used to probe the simulation
#           box
#
def surface_list(atoms, r_min, r_max, n_tries, cluster=False):
    if not cluster:
        raise ValueError("Currently the surface_list function can only be used for cluster systems")
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    is_surface = surface_tools.surface_atoms(pos, cell, r_min, r_max, n_tries, cluster)
    surf_list = []
    for i in range(0, len(atoms)):
        if is_surface[i] == 1:
            surf_list.append(i)
    return surf_list
###############################################################################################
