from ase.io import read
import numpy as np
from ase_tools_for_ovito import Pipeline


# Example plotting a GST structure with default parameters
#
#  - For rendering atoms and bonds:
#    - If radii are not given, the vdW radii from ASE times 0.3 will be used
#    - If bond cutoffs are not given, ASE's covalent radii times 2.5 will be used
#    - If colors are not given, the JMOL defaults will be used
#  - For the scene options:
#    - If no scale is given, it defaults to volume**(1/3) times 3 (the bigger the scale the smaller the atoms appear)
#    - If no direction is given, [3,2,-1] will be used
#    - If no shift is given, [0,0,0] will be used
#  - For exporting the file:
#    - The file name must always be specified by the user
#    - If alpha is not specified, it's set to True by default
#    - If no size is specified, the default is 512x512 pixels


# Create an ASE Atoms object
atoms = read("GST.xyz")

# Create a Pipeline object from atoms with all options set to defaults
pipeline = Pipeline(atoms)

# Set the visual elements with defaults
pipeline.set_visuals()

# Export to a PNG file
pipeline.render("defaults.png")




# Example modifying the various options in a "destructive way"
pipeline = Pipeline(atoms,
                    radii={"Ge": 0.1, "Sb": 0.2, "Te": 0.3}, # exaggeratedly small radii to check the effect
                    colors={"Ge": [1,0,0], "Sb": [0,1,0], "Te": [0,0,1]}, # RGB colors, Ge = red, Sb = green, Te = blue
                    cutoffs={"Ge": 0., "Sb": 3., "Te": 3.}) # to draw Sb-Te bonds but not Ge-* bonds

pipeline.set_visuals(direction=[1,1,1], # look at it from the corner of the unit cell
                     scale=1., # very close look at the system
                     shift=[0,0,0]) # no shift of the position

pipeline.render("custom1.png",
                size=(1024,512), # we make a wider aspect-ratio figure
                alpha=True) # keep alpha, why not



# Just change the colors
pipeline = Pipeline(atoms, colors={"Ge": [1,0,0], "Sb": [0,1,0], "Te": [0,0,1]})
pipeline.set_visuals()
pipeline.render("custom2.png")



# Just change the atom sizes
pipeline = Pipeline(atoms, radii={"Ge": 0.1, "Sb": 0.2, "Te": 0.3})
pipeline.set_visuals()
pipeline.render("custom3.png")



# Provide colors with an array (random here, but you might have actual useful info to display!)
pipeline = Pipeline(atoms, colors=np.random.sample([len(atoms), 3]))
pipeline.set_visuals()
pipeline.render("custom4.png")
