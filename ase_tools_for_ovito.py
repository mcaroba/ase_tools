# ase_tools is copyright (c) 2022-2025 of Miguel Caro and Rina Ibragimova
# See License.md for more information

import numpy as np
from ase.io import read
from ase.data import covalent_radii, vdw_radii, atomic_numbers
from ovito.io.ase import ase_to_ovito
from ovito.pipeline import StaticSource
from ovito.pipeline import Pipeline as OvitoPipeline
from ovito.data import DataCollection

###############################################################################################
#
class Pipeline:
    def __init__(self, atoms, radii=None):
        self.atoms = atoms # This is an ASE Atoms object
        # Ovito has trouble with a bunch of array properties assigned to Atoms objects. The fix_atoms
        # property assigned by TurboGAP needs to be removed!
        del self.atoms.arrays["fix_atoms"]
        #
        atoms_data = ase_to_ovito(self.atoms)
        pipeline = OvitoPipeline(source = StaticSource(data = atoms_data)) # This is an Ovito Pipeline
        data = pipeline.compute()
        # Build a dictionary of chemical elements
        el_num = {}
        for at in self.atoms:
            if at.symbol not in el_num:
                el_num[at.symbol] = None
        types = data.particles_.particle_types_
        for n in range(1, len(el_num)+1):
            el = types.type_by_id_(n).name
            el_num[el] = n
        self.el_num = el_num
        # Set atomic radii
        if radii is not None:
            self.radii = radii
        else:
            radii = {}
            for el in el_num:
                numb = atomic_numbers[el]
                radii[el] = vdw_radii[numb] * 0.3
            self.radii = radii

        # Modifier functions
        def set_particle_size(frame: int, data: DataCollection):
            for el in self.el_num:
                n = self.el_num[el]
                data.particles_.particle_types_.type_by_id_(n).radius = self.radii[el]
        pipeline.modifiers.append(set_particle_size)

        # Return Ovito pipeline
        self.pipeline = pipeline
