"""
soap_turbo.descriptor
=====================

High-level Python interface to the soap_turbo Fortran library.

Usage
-----
    from soap_turbo.descriptor import Descriptor
    import ase.io

    desc_H = Descriptor(
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
        compress_soap=True,
        compress_mode="trivial"
    )

    atoms = ase.io.read("structure.xyz")
    soap            = desc_H.calc(atoms)
    soap, soap_der  = desc_H.calc(atoms, derivatives=True)
"""

import numpy as np
from ase.neighborlist import neighbor_list

try:
    # Import the patched f90wrap wrapper directly rather than the package,
    # to avoid a circular import (soap_turbo/__init__.py imports this file).
    from soap_turbo import soap_turbo_ext as _soap_turbo_mod
    _get_soap             = _soap_turbo_mod.soap_turbo_desc.get_soap
    _get_compress_indices = _soap_turbo_mod.soap_turbo_compress_module.get_compress_indices
except ImportError:
    raise ImportError(
        "soap_turbo extension module not found. "
        "Build and install it with ./build.sh before importing this module."
    )


_Z_TO_SYM = {
    1: "H",  2: "He", 3: "Li", 4: "Be", 5: "B",  6: "C",  7: "N",  8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",  16: "S",
    17: "Cl", 18: "Ar", 19: "K",  20: "Ca", 22: "Ti", 24: "Cr", 25: "Mn",
    26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 32: "Ge", 47: "Ag",
    48: "Cd", 50: "Sn", 74: "W",  78: "Pt", 79: "Au", 82: "Pb",
}


class Descriptor:
    """
    A soap_turbo descriptor object.

    Parameters
    ----------
    n_max : int
        Number of radial basis functions per species.
    l_max : int
        Maximum angular momentum quantum number.
    alpha_max : list of int, length n_species
        Number of radial basis functions per species.
    atom_sigma_r : list of float, length n_species
        Radial width of atomic Gaussians per species.
    atom_sigma_t : list of float, length n_species
        Tangential width of atomic Gaussians per species.
    atom_sigma_r_scaling : list of float, length n_species
        Radial Gaussian width scaling per species.
    atom_sigma_t_scaling : list of float, length n_species
        Tangential Gaussian width scaling per species.
    rcut_hard : float
        Hard cutoff radius (Angstrom).
    rcut_soft : float
        Soft cutoff radius (Angstrom). Smooth decay begins here.
    basis : str
        Radial basis type: "poly3", "poly3gauss", or "poly3_tabulated".
    scaling_mode : str
        Amplitude scaling mode: "polynomial", "exponential", or "trivial".
    amplitude_scaling : list of float, length n_species
        Per-species amplitude scaling exponents.
    species_Z : list of int
        Atomic numbers of all species in the descriptor, in order.
        Order determines species indices (1-based in Fortran).
    central_index : int
        1-based index into species_Z for the central atom species.
    radial_enhancement : int
        Radial enhancement factor (0 or 1).
    central_weight : list of float, length n_species
        Weight of the central atom contribution per species.
    nf : list of float, length n_species, optional
        Normalisation factors. Defaults to 1.0 for all species.
    global_scaling : list of float, length n_species, optional
        Global scaling. Defaults to 1.0 for all species.
    compress_mode : str, optional
        Compression mode. Default "trivial" (no compression).
    do_timing : bool, optional
        Print Fortran-side timing. Default False.
    """

    def __init__(
        self,
        n_max,
        l_max,
        alpha_max,
        atom_sigma_r,
        atom_sigma_t,
        atom_sigma_r_scaling,
        atom_sigma_t_scaling,
        rcut_hard,
        rcut_soft,
        basis,
        scaling_mode,
        amplitude_scaling,
        species_Z,
        central_index,
        radial_enhancement,
        central_weight,
        nf=None,
        global_scaling=None,
        compress_soap=True,
        compress_mode="trivial",
        do_timing=False,
    ):
        self.n_max            = int(n_max)
        self.l_max            = int(l_max)
        self.n_species        = len(species_Z)
        self.species_Z        = list(species_Z)
        self.central_index    = int(central_index)   # 1-based
        self.rcut_hard        = float(rcut_hard)
        self.rcut_soft        = float(rcut_soft)
        self.basis            = str(basis)
        self.scaling_mode     = str(scaling_mode)
        self.compress_soap    = bool(compress_soap)
        self.compress_mode    = str(compress_mode)
        self.radial_enhancement = int(radial_enhancement)
        self.do_timing        = bool(do_timing)

        ns = self.n_species

        def _arr(x, name, dtype=float):
            a = np.asarray(x, dtype=dtype)
            if a.shape != (ns,):
                raise ValueError(
                    f"{name} must have length n_species={ns}, got shape {a.shape}"
                )
            return np.asfortranarray(a)

        self.alpha_max            = _arr(alpha_max,             "alpha_max", int)
        self.atom_sigma_r         = _arr(atom_sigma_r,          "atom_sigma_r")
        self.atom_sigma_t         = _arr(atom_sigma_t,          "atom_sigma_t")
        self.atom_sigma_r_scaling = _arr(atom_sigma_r_scaling,  "atom_sigma_r_scaling")
        self.atom_sigma_t_scaling = _arr(atom_sigma_t_scaling,  "atom_sigma_t_scaling")
        self.amplitude_scaling    = _arr(amplitude_scaling,     "amplitude_scaling")
        self.central_weight       = _arr(central_weight,        "central_weight")
        self.nf             = _arr(nf if nf is not None else [1.0]*ns, "nf")
        self.global_scaling = _arr(global_scaling if global_scaling is not None else [1.0]*ns, 
                                    "global_scaling")
        self.rcut_hard_arr  = np.asfortranarray(np.full(ns, self.rcut_hard))
        self.rcut_soft_arr  = np.asfortranarray(np.full(ns, self.rcut_soft))

        self._setup_compression()
        #self._compression_indices()


    # ------------------------------------------------------------------ #
    #  Compression                                                         #
    # ------------------------------------------------------------------ #

    # Get compression indices
    # CURRENTLY NOT IN USE
    def _compression_indices(self):
        # Set up
        n_max = 0
        for i in range(self.n_species):
            n_max += self.alpha_max[i]
        #n_max = np.sum(self.alpha_max)
        #max_size = (n_max * (n_max + 1) // 2 * (self.l_max + 1)) ** 2 + 1
        self._n_soap_uncompressed = n_max * (n_max + 1) // 2 * (self.l_max + 1)
        P_i = np.zeros(self._n_soap_uncompressed, dtype=np.int32, order='F')
        P_j = np.zeros(self._n_soap_uncompressed, dtype=np.int32, order='F')
        P_el = np.zeros(self._n_soap_uncompressed, dtype=np.float64, order='F')
        if self.compress_soap:
            # Pythonized trivial compression algorithm
            if self.compress_mode == "trivial" or self.compress_mode == "Trivial":
                pivot = np.zeros(self.n_species, dtype=int)
                pivot[0] = 1
                for i in range(len(self.alpha_max)-1):
                    pivot[i+1] = pivot[i] + self.alpha_max[i]
                
                counter = 0
                k = 1
                for n in range(1, n_max + 1):
                    for m in range(n, n_max + 1):
                        for l in range(0, self.l_max + 1):
                            if np.any(n == pivot) or np.any(m == pivot):
                                counter += 1
                                P_i[counter - 1] = counter
                                P_j[counter - 1] = k
                                P_el[counter - 1] = 1.0
                            k += 1
                
                self._compress_P_nonzero = counter
                self._n_soap = counter
                self._compress_P_i = np.asfortranarray(P_i[:self._n_soap].copy())
                self._compress_P_j = np.asfortranarray(P_j[:self._n_soap].copy())
                self._compress_P_el = np.asfortranarray(P_el[:self._n_soap].copy())

        else:
            self._n_soap = self._n_soap_uncompressed
            self._compress_P_i = np.asfortranarray(P_i.copy())
            self._compress_P_j = np.asfortranarray(P_i.copy())
            self._compress_P_el = np.asfortranarray(P_i.copy())
            self._compress_P_nonzero = 0
            
        
    # Pass copies of every array argument to _get_compress_indices.
    def _setup_compression(self):
        """
        Call get_compress_indices to determine output SOAP dimension and
        fill compression index arrays.

        Available compression methods:
            None
            Trivial
            0_0
            0_1
            0_2
            1_0
            1_1
            1_2
            2_0
            2_1
            2_2
        """

        n_max       = int(np.max(self.alpha_max))
        max_size    = (n_max * (n_max + 1) // 2 * (self.l_max + 1)) ** 2 + 1

        P_i   = np.zeros(max_size, dtype=np.int32,   order='F')
        P_j   = np.zeros(max_size, dtype=np.int32,   order='F')
        P_el  = np.zeros(max_size, dtype=np.float64, order='F')
        dim  = np.array([0], dtype=np.int32)
        P_nz = np.array([0], dtype=np.int32)

        if self.compress_soap:
            # First call: determine dimension only
            _get_compress_indices(
                self.compress_mode,
                self.alpha_max,
                self.l_max,
                dim,
                P_nz,
                P_i,
                P_j,
                P_el,
                "get_dim",
            )

            # Second call: fill index arrays
            _get_compress_indices(
                self.compress_mode,
                self.alpha_max,
                self.l_max,
                dim,
                P_nz,
                P_i,
                P_j,
                P_el,
                "set_indices",
            )

        n = int(P_nz[0])
        self._n_soap = int(dim[0])
        self._compress_P_nonzero = int(P_nz[0])
        self._compress_P_i   = np.asfortranarray(P_i[:n].copy())
        self._compress_P_j   = np.asfortranarray(P_j[:n].copy())
        self._compress_P_el  = np.asfortranarray(P_el[:n].copy())

        # Uncompressed dimension — n_max = SUM(alpha_max), not max.
        self._n_soap_uncompressed = n_max * (n_max + 1) // 2 * (self.l_max + 1)
        if not self.compress_soap:
            self._n_soap = self._n_soap_uncompressed

    # ------------------------------------------------------------------ #
    #  Neighbour list                                                      #
    # ------------------------------------------------------------------ #

    def _build_neighbour_arrays(self, atoms):
        """
        Build all per-pair arrays required by get_soap.

        Array shapes (Fortran column-major, so shape here is Python row-major):
          rjs          : (n_atom_pairs,)
          thetas       : (n_atom_pairs,)
          phis         : (n_atom_pairs,)
          n_neigh      : (n_sites,) integer
          mask         : (n_atom_pairs, n_species)  logical, Fortran order
          species      : (max_mult, n_sites) integer, Fortran order
          species_multiplicity : (n_sites,)  integer
        """
        central_Z    = self.species_Z[self.central_index - 1]
        all_Z        = atoms.get_atomic_numbers()
        central_sites = [i for i, Z in enumerate(all_Z) if Z == central_Z]

        if not central_sites:
            raise ValueError(
                f"No atoms of species Z={central_Z} "
                f"(central_index={self.central_index}) found."
            )

        n_sites = len(central_sites)

        # ASE neighbour list for ALL atom pairs within rcut_hard
        i_idx, j_idx, d_vec, dist = neighbor_list("ijDd", atoms, self.rcut_hard)

        # Build per-site data in central-site order
        rjs_list, thetas_list, phis_list = [], [], []
        n_neigh = np.zeros(n_sites, dtype=np.int32)

        for site_pos, site_i in enumerate(central_sites):
            sel   = i_idx == site_i
            vecs  = d_vec[sel]
            dists = dist[sel]
            n_neigh[site_pos] = len(dists)

            r     = dists
            x, y, z = vecs[:, 0], vecs[:, 1], vecs[:, 2]
            theta = np.arccos(np.clip(z / np.where(r > 0, r, 1.0), -1.0, 1.0))
            phi   = np.arctan2(y, x)

            rjs_list.append(dists)
            thetas_list.append(theta)
            phis_list.append(phi)

        # Fortran expects the central atom itself as neighbour j=1 for each site
        # (at distance rj=0). It is skipped in the radial expansion unless
        # do_central=.true., but it must be present in the list so that all
        # subsequent neighbour indices are correct. Prepend a self-pair to each
        # site's neighbour list and increment n_neigh by 1 per site.
        for site_pos in range(n_sites):
            rjs_list[site_pos]    = np.concatenate([[0.0],    rjs_list[site_pos]])
            thetas_list[site_pos] = np.concatenate([[0.0],    thetas_list[site_pos]])
            phis_list[site_pos]   = np.concatenate([[0.0],    phis_list[site_pos]])
            n_neigh[site_pos]    += 1

        rjs    = np.concatenate(rjs_list)    if rjs_list    else np.empty(0)
        thetas = np.concatenate(thetas_list) if thetas_list else np.empty(0)
        phis   = np.concatenate(phis_list)   if phis_list   else np.empty(0)
        n_atom_pairs = int(rjs.size)

        # species(max_mult, n_sites) — Fortran order
        # For now every site has exactly one species (multiplicity 1).
        species = np.zeros((1, n_sites), dtype=np.int32, order='F')
        species_multiplicity = np.ones(n_sites, dtype=np.int32)
        for site_pos, site_i in enumerate(central_sites):
            Z = all_Z[site_i]
            species[0, site_pos] = self.species_Z.index(Z) + 1   # 1-based

        # mask(n_atom_pairs, n_species) — Fortran order, int32 (not bool).
        # The first pair for each site is the self-pair; set its mask to the
        # central atom's own species. Then set masks for external neighbours.
        mask = np.zeros((n_atom_pairs, self.n_species), dtype=np.int32, order='F')
        pair_idx = 0
        central_sp_idx = self.species_Z.index(self.species_Z[self.central_index - 1])
        for site_i in central_sites:
            # Self-pair: species of the central atom itself
            mask[pair_idx, central_sp_idx] = 1
            pair_idx += 1
            sel      = i_idx == site_i
            neigh_js = j_idx[sel]
            for jj in neigh_js:
                Z_j = all_Z[jj]
                if Z_j in self.species_Z:
                    sp_idx = self.species_Z.index(Z_j)   # 0-based
                    mask[pair_idx, sp_idx] = 1
                pair_idx += 1

        return (
            np.asfortranarray(rjs.astype(np.float64)),
            np.asfortranarray(thetas.astype(np.float64)),
            np.asfortranarray(phis.astype(np.float64)),
            np.asfortranarray(n_neigh),
            mask,
            species,
            np.asfortranarray(species_multiplicity),
            n_sites,
            n_atom_pairs,
        )

    # ------------------------------------------------------------------ #
    #  Public API                                                          #
    # ------------------------------------------------------------------ #

    def calc(self, atoms, derivatives=False):
        """
        Compute SOAP descriptors for all central-species atoms in `atoms`.

        Parameters
        ----------
        atoms : ase.Atoms
            The structure. Set pbc=True and a cell for periodic structures.
        derivatives : bool, optional
            If True, also return Cartesian derivatives. Default False.

        Returns
        -------
        soap : np.ndarray, shape (n_central_atoms, n_soap)
            SOAP descriptor vectors, one row per central atom.
        soap_cart_der : np.ndarray, shape (3, n_soap, n_atom_pairs)
            Cartesian derivatives. Only returned when derivatives=True.
        """
        (rjs, thetas, phis, n_neigh, mask, species,
         species_multiplicity, n_sites, n_atom_pairs) = \
            self._build_neighbour_arrays(atoms)

        # Fortran array layouts (column-major):
        #   soap          : (n_soap, n_sites)
        #   soap_cart_der : (3, n_soap, n_atom_pairs)
        soap_buf     = np.zeros((self._n_soap, n_sites),
                                dtype=np.float64, order='F')
        soap_der_buf = np.zeros((3, self._n_soap, n_atom_pairs),
                                dtype=np.float64, order='F')

        _get_soap(
            n_sites              = n_sites,
            n_neigh              = n_neigh,
            n_species            = self.n_species,
            species              = species,
            species_multiplicity = species_multiplicity,
            n_atom_pairs         = n_atom_pairs,
            mask                 = mask,
            rjs                  = rjs,
            thetas               = thetas,
            phis                 = phis,
            alpha_max            = self.alpha_max,
            l_max                = self.l_max,
            rcut_hard            = self.rcut_hard_arr,
            rcut_soft            = self.rcut_soft_arr,
            nf                   = self.nf,
            global_scaling       = self.global_scaling,
            atom_sigma_r         = self.atom_sigma_r,
            atom_sigma_r_scaling = self.atom_sigma_r_scaling,
            atom_sigma_t         = self.atom_sigma_t,
            atom_sigma_t_scaling = self.atom_sigma_t_scaling,
            amplitude_scaling    = self.amplitude_scaling,
            radial_enhancement   = self.radial_enhancement,
            central_weight       = self.central_weight,
            basis                = self.basis,
            scaling_mode         = self.scaling_mode,
            do_timing            = self.do_timing,
            do_derivatives       = bool(derivatives),
            compress_soap        = self.compress_soap,
            compress_p_nonzero   = self._compress_P_nonzero,
            compress_p_i         = self._compress_P_i,
            compress_p_j         = self._compress_P_j,
            compress_p_el        = self._compress_P_el,
            soap                 = soap_buf,
            soap_cart_der        = soap_der_buf,
        )

        # Copy results to new Python-owned arrays.
        if derivatives:
            return soap_buf.T, soap_der_buf.T

        return soap_buf.T

    # ------------------------------------------------------------------ #
    #  Repr                                                                #
    # ------------------------------------------------------------------ #

    def __repr__(self):
        central_sym  = _Z_TO_SYM.get(self.species_Z[self.central_index - 1],
                                      str(self.species_Z[self.central_index - 1]))
        species_syms = [_Z_TO_SYM.get(Z, str(Z)) for Z in self.species_Z]
        return (
            f"Descriptor("
            f"central={central_sym}, "
            f"species={species_syms}, "
            f"n_max={self.n_max}, l_max={self.l_max}, "
            f"rcut_hard={self.rcut_hard}, rcut_soft={self.rcut_soft}, "
            f"basis={self.basis!r}, "
            f"n_soap={self._n_soap}"
            f")"
        )
