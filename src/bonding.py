#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = 'Chao Yang'
__version__=	'1.0'



"""
Define the bonding related functions.

Electronic negativity strategy:
    Refer to Matterviz developed by Janosh.
    * Electronegativity-based bonding with chemical preferences.
    * This algorithm considers electronegativity differences between atoms, metal/nonmetal
    * properties, and distance to determine bond strength. Bonds are only created if the
    * computed strength exceeds the strength_threshold parameter (default: 0.3).

"""

from pymatgen.core import Site, Structure


def build_spatial_grid(sites: list[Site], cell_size: float) -> dict[tuple[int, int, int], list[int]]:
    """
    Build a spatial grid for efficient neighbor searching.

    Parameters:
    -----------
    sites: list of Site objects
        The atomic sites in the structure.
    cell_size: float
        The size of each grid cell.
    Returns:
    --------
    grid: dict [tuple[int, int, int], list[int]]
        A dictionary mapping grid cell indices to lists of site indices.
    """
    grid = {}
    for i, site in enumerate(sites):
        x, y, z = site.coords
        ix = int(x // cell_size)
        iy = int(y // cell_size)
        iz = int(z // cell_size)
        key = (ix, iy, iz)
        if key not in grid:
            grid[key] = []
        grid[key].append(i)
    return grid

def electroneg_ratio(
        structure: Structure,
        electroneg_threshold: float = 1.7, # Maximum electronegativity difference for bonding
        max_distance_ratio: float = 2.0, # Factor for the sum of covalent radii
        min_bond_distance: float = 0.4, # Minimum bond distance in Angstrom
        metal_metal_penalty: float = 0.5, # Strength penalty for metal-metal bonding
        metal_nonmetal_bonus: float = 1.5, # Strength bonus for metal-nonmetal bonding
        similarity_penalty: float = 1.2, # Bonus for similar electronegativity
        same_element_penalty: float = 0.5, # Penalty for same element bonding
        strengthen_bonus: float = 0.3 # Minimum bonding strength to include in results
    ):
    """
    Calculate the bonding status based on the electronegativity ratio.
    """

    sites = structure.sites
    num_sites = len(sites)


def main():
    Structure.from_file("working/surface_model/surface/miller_010_layers2_minlat10.0_symTrue/slab_1.cif")


if __name__ == '__main__':
    main()
    




