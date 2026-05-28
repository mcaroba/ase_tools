"""
Test / example script for soap_turbo.descriptor.

Computes SOAP descriptors for H and C atoms in a small hydrocarbon structure,
mirroring the quippy usage pattern.
"""

import numpy as np
import ase
from ase import Atoms
from ase.io import write

# ── Test import ────────────────────────────────────────────────────────────
import soap_turbo
print("Loaded:", soap_turbo)
print("Contains:", dir(soap_turbo))

# ── Build a small test structure: methane CH4 ──────────────────────────────
d = 1.089   # C-H bond length in Angstrom
atoms = Atoms(
    symbols=["C", "H", "H", "H", "H"],
    positions=[
        [0.000,  0.000,  0.000],
        [d,      d,      d    ],
        [-d,    -d,      d    ],
        [-d,     d,     -d    ],
        [d,     -d,     -d    ],
    ],
    cell=[30.0, 30.0, 30.0],
    pbc=True
)

# ── Define descriptors for each species ────────────────────────────────────
common = dict(
    n_max=8, l_max=8,
    alpha_max=[8, 8],
    atom_sigma_r=[0.2, 0.2],
    atom_sigma_t=[0.2, 0.2],
    atom_sigma_r_scaling=[0.1, 0.1],
    atom_sigma_t_scaling=[0.1, 0.1],
    rcut_hard=5.0,
    rcut_soft=4.5,
    basis="poly3gauss",
    scaling_mode="polynomial",
    amplitude_scaling=[2.0, 2.0],
    species_Z=[1, 6],       # H=1, C=6
    radial_enhancement=1,
    central_weight=[1.0, 1.0],
    compress_soap=True,
    compress_mode="Trivial"
)

desc_H = soap_turbo.Descriptor(**common, central_index=1)   # descriptors centred on H
desc_C = soap_turbo.Descriptor(**common, central_index=2)   # descriptors centred on C

print("Descriptor for H:", desc_H)
print("Descriptor for C:", desc_C)
print()

# ── Compute descriptors ────────────────────────────────────────────────────
soap_H = desc_H.calc(atoms)
soap_C = desc_C.calc(atoms)

print(f"SOAP for H atoms: shape = {soap_H.shape}")   # (n_H_atoms, n_soap)
print(f"SOAP for C atoms: shape = {soap_C.shape}")   # (n_C_atoms, n_soap)
print()
print("Testing descriptor lengths:")
print("H descriptor lengths:", np.linalg.norm(soap_H, axis=1))
print("C descriptor lengths:", np.linalg.norm(soap_C, axis=1))

# ── With derivatives ───────────────────────────────────────────────────────
#soap_H, dsoap_H = desc_H.calc(atoms, derivatives=True)
#print(f"SOAP derivatives for H: shape = {dsoap_H.shape}")
