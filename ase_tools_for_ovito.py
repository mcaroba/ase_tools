# ase_tools is copyright (c) 2022-2025 of Miguel Caro and Rina Ibragimova
# See License.md for more information

import numpy as np
from ase.io import read
from ase.data import covalent_radii, vdw_radii, atomic_numbers
from ase.data.colors import jmol_colors, cpk_colors
from ovito.io.ase import ase_to_ovito
from ovito.pipeline import StaticSource
from ovito.pipeline import Pipeline as OvitoPipeline
from ovito.data import DataCollection
from ovito.modifiers import CreateBondsModifier
from ovito.vis import Viewport, TachyonRenderer

###############################################################################################
#
class Pipeline:
    def __init__(self, atoms, radii=None, cutoffs=None, colors=None):
        self.atoms = atoms # This is an ASE Atoms object

        # Ovito has trouble with a bunch of array properties assigned to Atoms objects. The fix_atoms
        # property assigned by TurboGAP needs to be removed!
        del self.atoms.arrays["fix_atoms"]

        # Create the data objects
        atoms_data = ase_to_ovito(self.atoms)
        pipeline = OvitoPipeline(source = StaticSource(data = atoms_data)) # This is an Ovito Pipeline
        data = pipeline.compute()

        # Build various system properties

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

        # Set bonding cutoffs
        if cutoffs is not None:
            self.cutoffs = cutoffs
        else:
            cutoffs = {}
            for el in el_num:
                numb = atomic_numbers[el]
                cutoffs[el] = covalent_radii[numb] * 2.5
            self.cutoffs = cutoffs

        # Set various self properties
        if type(colors) is str:
            self.colors = colors.lower()
        else:
            self.colors = colors
        self.mincolor = np.array([0., 0., 1.])
        self.maxcolor = np.array([1., 0., 0.])

        # Store Ovito pipeline in self to be accessible by modifiers
        self.pipeline = pipeline

        # Append modifiers
        self.set_particle_size()
        self.create_bonds()
        self.set_particle_color()

    # >>>>>>> Modifier functions
    # TODO: allow for modifier-specific parameter modification (let these functions take
    #       parameters as input rather than retrieve all of them from the class
    #       instantiation)

    # Set the size of the atoms
    def set_particle_size(self):
        def modifier(frame: int, data: DataCollection):
            for el in self.el_num:
                n = self.el_num[el]
                data.particles_.particle_types_.type_by_id_(n).radius = self.radii[el]
        self.pipeline.modifiers.append(modifier)

    # Set the cutoff to draw the bonds
    def create_bonds(self):
        modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise) # This uses Ovito's own modifier
        for el in self.el_num:
            n = self.el_num[el]
            for el2 in self.el_num:
                n2 = self.el_num[el2]
                if n2 >= n:
                    modifier.set_pairwise_cutoff(n, n2, (self.cutoffs[el]+self.cutoffs[el2])/2.)
        self.pipeline.modifiers.append(modifier)

    # Color particles
    def set_particle_color(self):
        if self.colors == None or self.colors == "jmol": # JMOL color scheme (ASE's default)
            colors = {}
            for el in self.el_num:
                colors[el] = jmol_colors[atomic_numbers[el]]
            self.colors = colors
        elif self.colors == "cpk": # CPK color scheme
            colors = {}
            for el in self.el_num:
                colors[el] = cpk_colors[atomic_numbers[el]]
            self.colors = colors
        if type(self.colors) is dict: # If self.colors is a dictionary, we create the color array from it
            colors = np.zeros([len(self.atoms), 3])
            for i in range(0, len(self.atoms)):
                colors[i] = self.colors[self.atoms[i].symbol]
            self.colors = colors
        elif type(self.colors) is str: # If self.colors is a string, we assume it refers to an array property
            colors = np.zeros([len(self.atoms), 3])
            array = self.atoms.get_array(self.colors)
            amin = np.min(array)
            amax = np.max(array)
            for i in range(0, len(self.atoms)):
                if amax - amin > 0.:
                    colors[i] = self.mincolor + (self.maxcolor - self.mincolor) / (amax - amin) * (array[i] - amin)
                else:
                    colors[i] = self.mincolor
            self.colors = colors
        elif len(self.colors) == len(self.atoms): # Otherwise, we assume the user is providing a list or array
            pass
        else:
            raise("Bad colors option!")

        def modifier(frame: int, data: DataCollection):
            cols = data.particles_.create_property('Color')
            for i in range(0, len(cols)):
                cols[i] = self.colors[i]
        self.pipeline.modifiers.append(modifier)

    # End of modifier functions <<<<<<<

    # Visualization and rendering options
    def set_visuals(self, renderer="tachyon",
                    shadows=False, direct_light_intensity=1.1,
                    direction=[1,1,1], scale=1., shift=[0,0,0]):
        if renderer == "tachyon":
            self.renderer = TachyonRenderer(shadows=shadows, direct_light_intensity=direct_light_intensity)
        else:
            raise("%s not supported as renderer option")
        self.vp = Viewport()
        self.vp.type = Viewport.Type.Perspective
        self.vp.camera_dir = direction
        self.vp.camera_pos = shift + self.atoms.get_center_of_mass() - direction / np.dot(direction,direction)**0.5 * scale

    # Trigger rendering
    def render(self, filename, size=(512,512), alpha=True):
        self.pipeline.add_to_scene()
        self.vp.render_image(filename = filename, size=size, alpha=alpha, renderer=self.renderer)
        self.pipeline.remove_from_scene()
