from ase.visualize import view
from ase import Atoms
from ase.io import read
atoms = read("graphene_C60.xyz")
from cluster import cluster_module as cl
import numpy as np
pos = atoms.get_positions()
cell = atoms.get_cell()
lx = cell[0][0]
ly = cell[1][1]
lz = cell[2][2]
list = cl.cluster_atoms(np.transpose(pos), lx, ly, lz, 1.8)


pos_new = []
for i in range(0, len(list)):
    if list[i] == 1:
        pos_new.append(pos[i])

mol1 = Atoms("C%i" % len(pos_new), cell = cell, positions=pos_new, pbc=True)

pos_new = []
for i in range(0, len(list)):
    if list[i] == 2:
        pos_new.append(pos[i])

mol2 = Atoms("C%i" % len(pos_new), cell = cell, positions=pos_new, pbc=True)

view(mol1)
view(mol2)
