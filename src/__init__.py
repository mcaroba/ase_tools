"""
soap_turbo — Python interface to the soap_turbo Fortran library.

Usage
-----
    from soap_turbo import Descriptor
    import ase.io

    desc = Descriptor(
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
        species_Z=[1, 6],
        central_index=1,
        radial_enhancement=1,
        central_weight=[1.0, 1.0],
    )

    atoms = ase.io.read("structure.xyz")
    soap = desc.calc(atoms)
"""

from .descriptor import Descriptor

__all__ = ["Descriptor"]
