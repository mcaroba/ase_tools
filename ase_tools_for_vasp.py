# ase_tools is copyright (c) 2022-2023 of Miguel Caro and Rina Ibragimova
# See License.md for more information

import numpy as np
from ase.io import read

###############################################################################################
#
# This function makes a VASP KPOINTS files compatible with the chosen combination of KSPACING
# and KGAMMA parameters. You can also apply a shift to the k-mesh origin. The main
# functionality is that it can detect the presence of vacuum in the simulation cell, such
# that the sampling along the vacuum direction(s) will be one k point only. This can help
# to save CPU time in situations where a high-throughput approach is being used and the user
# wants to avoid unnecessarily dense k-point sampling for surfaces, nanoparticles, etc.
#
# The vacuum_check variable enables detection of vacuum, and the vacuum_cutoff variable tells
# the code what gap in the atomic density corresponds to the presence of vacuum. E.g.,
# vacuum_cutoff = 5. means that whenever there is a 5 Angstrom gap with no atoms present
# along some spatial direction, the corresponding k-sampling will be set to 1 along that
# direction. When this leads to Gamma-point sampling (e.g., for nanoparticles and molecules),
# the user should further use this information to use the Gamma-point VASP binary, for extra
# CPU time savings.
#
# Currently, the code can only detect trivial vacuum gaps perpendicular to a given plane (where
# this plane is defined by two of the lattice vectors through their cross product). Complex
# "curved" gaps cannot be detected. This should not be an issue for most applications.
#
def make_kpoints_file(atoms, filename="KPOINTS", kspacing=0.5, kgamma=True,
                      shift=[0.,0.,0.], vacuum_check=True, vacuum_cutoff=5.):
    at = atoms.copy()
    at.wrap()
    cell = at.get_cell()
    pos = at.get_positions()
    a = np.zeros([3,3])
    a[0] = cell[0,:]; a[1] = cell[1,:]; a[2] = cell[2,:]
    u = np.zeros([3,3])
    u[0] = np.cross(a[1], a[2]); u[1] = np.cross(a[2], a[0]); u[2] = np.cross(a[0], a[1])
    vol = np.dot(a[0], u[0])
    # The convention is: ai dot bj = dij (no 2pi factor)
    b = np.zeros([3,3])
    b[0] = u[0]/vol; b[1] = u[1]/vol; b[2] = u[2]/vol
    # VASP's correspondence between KSPACING and automatic KPOINTS generation:
    nk = [0, 0, 0]
    for i in range(0, 3):
        nk[i] = round(max(1., 2.*np.pi/kspacing*np.dot(b[i],b[i])**0.5 + 0.5))
    if vacuum_check:
        for i in range(0, 3):
            more_pos = []
            for k in [0, 1]:
                for j in range(0, len(pos)):
                    more_pos.append(pos[j] + float(k)*a[i])
            proj = np.dot(more_pos, u[i]/np.dot(u[i],u[i])**0.5)
            proj = np.sort(proj)
            max_dist = 0.
            for j in range(0, len(proj)-1):
                if proj[j+1] - proj[j] > max_dist:
                    max_dist = proj[j+1] - proj[j]
            if max_dist > vacuum_cutoff:
                nk[i] = 1
    # Write out
    f = open(filename, "w")
    print("Generated with ase_tool's make_kpoints_file()", file=f)
    print("0", file=f)
    if kgamma:
        print("Gamma", file=f)
    else:
        print("Monkhorst-Pack", file=f)
    print(*nk, file=f)
    print(*shift, file=f)
    f.close()
###############################################################################################




###############################################################################################
#
# This function provides a more powerful wrapper for ASE's native OUTCAR reader.
#
def read_vasp(outcar, check_convergence=True, extract_hirshfeld=True):
#   Use ASE reader to initialize the Atoms() object
    atoms = read(outcar)
#   Regular Python's read to import the OUTCAR's contents
    f = open(outcar, "r")
    lines = f.readlines()
    f.close()
#   Check that the electronic SCF is converged
    if check_convergence:
        n = 0
        nelm = None
        toten_list = []
        for line in lines:
            if len(line.split()) > 0:
                if line.split()[0] == "NELM":
                    nelm = int(line.split()[2][0:-1])
            if any([word == "TOTEN" for word in line.split()]):
                n += 1
                toten_list.append(float(line.split()[4]))
        if not nelm:
            raise RuntimeError("Your input file does not appear to be a valid VASP OUTCAR file!")
        if len(toten_list) >= 3:
            ediff = np.abs(toten_list[-3]-toten_list[-2])
            atoms.info["scf_ediff"] = ediff
        if n >= nelm+1:
            raise RuntimeError("Your VASP SCF calculation appears to not be converged " + \
                               "(the allowed maximum of NELM iterations was reached)! " + \
                               "You can disregard this error (at your own risk!) by setting " + \
                               "check_convergence=False. If you choose to do so, the energy " + \
                               "difference between the last two SCF steps (%.5e eV in this case)" % (ediff) + \
                               "will be returned with tag scf_ediff in the Atoms object's info dictionary.")
#   If a TS or MBD calculation was done, get Hirshfeld parameters
    if extract_hirshfeld:
        ivdw = 0
        it_hirshfeld = False
        for line in lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "IVDW":
                if line.split()[2] == "2" or line.split()[2] == "20":
                    ivdw = 20
                elif line.split()[2] == "21":
                    ivdw = 21
                    it_hirshfeld = True
                elif line.split()[2] == "202":
                    ivdw = 202
                elif line.split()[2] == "212":
                    ivdw = 212
                    it_hirshfeld = True
                break
        if ivdw != 0:
            hirsh_v, hirsh_q = extract_hirshfeld_volumes(lines)
            atoms.set_array("hirshfeld_v", hirsh_v)
            if len(hirsh_q) >= 1:
                atoms.set_array("hirshfeld_q", hirsh_q[0])
            if len(hirsh_q) == 2:
                atoms.set_array("hirshfeld_q_it", hirsh_q[1])
    return atoms

def extract_hirshfeld_volumes(lines):
    v = []
    q = []
    q_it = []
    for i in range(0, len(lines)):
        if lines[i].split() == ["Parameters", "of", "vdW", "forcefield:"]:
            j = i + 4
            while True:
                if len(lines[j].split()) == 0:
                    break
                else:
                    v.append(float(lines[j].split()[4]))
                    j += 1
        if lines[i].split() == ["Hirshfeld", "charges:"]:
            j = i + 3
            while True:
                if len(lines[j].split()) == 0:
                    break
                else:
                    q.append(float(lines[j].split()[2]))
                    j += 1
        if lines[i].split() == ["Hirshfeld-I", "charges:"]:
            j = i + 3
            while True:
                if len(lines[j].split()) == 0:
                    break
                else:
                    q_it.append(float(lines[j].split()[2]))
                    j += 1
    if len(q_it) > 0:
        return np.array(v), [np.array(q), np.array(q_it)]
    else:
        return np.array(v), [np.array(q)]
###############################################################################################
