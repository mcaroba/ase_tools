from ase.io import read
from ase_tools import split_atoms

atoms = read("atoms.xyz")
db = split_atoms(atoms, {"H": 1.3, "C": 1.8})
