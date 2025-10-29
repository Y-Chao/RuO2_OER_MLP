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

import math as m
from typing import Union
import numpy as np
from ase.data import covalent_radii
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

def setup_spatial_grid(
        sites: list[Site],
        cutoff: float
    ) -> tuple[Union[dict[tuple[int, int, int], list[int]], None], Union[float, None]]:
    """
    Setup the spatial grid for neighbor searching.

    Parameters:
    -----------
    sites: list of Site objects
        The atomic sites in the structure.
    cutoff: float
        The maximum cutoff distance for bonding.
    Returns:
    --------
    spatial_grid: dict [tuple[int, int, int], list[int]]
        The spatial grid for neighbor searching.
    cell_size: float or None
        The cell size used for the grid, or None if grid is not used.
    """
    if len(sites) > 50:
        return None, None
    return build_spatial_grid(sites, cutoff), cutoff

def get_neighbors_from_grid(
        pos: tuple[float, float, float],
        grid: dict[tuple[int, int, int], list[int]],
        cell_size: float,
    ) -> list[int]:
    """
    Get neighboring site indices from the spatial grid.
    Parameters:
    -----------
    pos: tuple[float, float, float]
        The position to search around.
    grid: dict[tuple[int, int, int], list[int]]
        The spatial grid.
    cell_size: float
        The size of each grid cell.
    Returns:
    --------
    neighbors: list[int]
        List of neighboring site indices.
    """
    [cx, cy, cz] = [m.floor(pos[0] / cell_size), m.floor(pos[1] / cell_size), m.floor(pos[2] / cell_size)]
    neighbors = []
    # Iterate over neighboring cells
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                cell_key = (cx + dx, cy + dy, cz + dz)
                if cell_key in grid:
                    neighbors.extend(grid[cell_key])
    return neighbors
    

def get_neighboring_candidates(
        pos: tuple[float, float, float],
        sites: list[Site],
        spatial_grid: Union[dict[tuple[int, int, int], list[int]], None],
        cell_size: Union[float, None],
    ) -> list[int]:
    """
    Get neighboring candidate site indices from the spatial grid.
    """
    if spatial_grid is None or cell_size is None:
        return list(range(len(sites)))
    return get_neighbors_from_grid(pos, spatial_grid, cell_size)

def compute_bond_transform(
        pos1: tuple[float, float, float],
        pos2: tuple[float, float, float]
    ) -> list[list[float]]:
    """
    Compute the transformation matrix for the bond between two positions.
    """
    dx, dy, dz = pos2[0] - pos1[0], pos2[1] - pos1[1], pos2[2] - pos1[2]
    height = m.hypot(dx, dy, dz) # Euclidean distance

    if height < 1e-10:
        return np.float32(np.eye(4))
    
    dir_x, dir_y, dir_z = dx / height, dy / height, dz / height
    [m00, m01, m02, m10, m11, m12, m20, m21, m22] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    # Special case: bond pointing straight up (+Y)
    if abs(dir_y - 1.0) < 1e-10:
        [m00, m01, m02, m10, m11, m12, m20, m21, m22] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    elif abs(dir_y + 1.0) < 1e-10:
        [m00, m01, m02, m10, m11, m12, m20, m21, m22] = [1, 0, 0, 0, -1, 0, 0, 0, 1]
    else:
        # General case: construct orthonormal basis (right, dir, up)
        # Right vector: perpendicular to dir in XZ plane
        [rx, rz] = [-dir_z, dir_x]
        right_length = m.hypot(rx, rz)
        [right_x, right_z] = [rx / right_length, rz / right_length]
        # Up vector: cross product of dir and right
        [up_x, up_y, up_z] = [
            dir_y * right_z,
            dir_z * right_x - dir_x * right_z,
            -dir_y * right_x
        ]
        [m00, m01, m02, m10, m11, m12, m20, m21, m22] = [
            right_x, dir_x, up_x,
            0, dir_y, up_y,
            right_z, dir_z, up_z
        ]
    
    # Position at midpoint between pos1 and pos2
    [px, py, pz] = [0.5 * (pos1[0] + pos2[0]),
                      0.5 * (pos1[1] + pos2[1]),
                      0.5 * (pos1[2] + pos2[2])]
    return np.float32([
        [m00, m10, m20, 0],
        [m01 * height, m11 * height, m21 * height, 0],
        [m02, m12, m22, 0],
        [px, py, pz, 1]
    ])
    

def electroneg_ratio(
        structure: Structure,
        electroneg_threshold: float = 1.7, # Maximum electronegativity difference for bonding
        max_distance_ratio: float = 2.0, # Factor for the sum of covalent radii
        min_bond_distance: float = 0.4, # Minimum bond distance in Angstrom
        metal_metal_penalty: float = 0.5, # Strength penalty for metal-metal bonding
        metal_nonmetal_bonus: float = 1.5, # Strength bonus for metal-nonmetal bonding
        similarity_x_bonus: float = 1.2, # Bonus for similar electronegativity
        same_element_penalty: float = 0.5, # Penalty for same element bonding
        strength_threshold: float = 0.3 # Minimum bonding strength to include in results
    ):
    """
    Calculate the bonding status based on the electronegativity ratio.
    """

    min_dist_sq = min_bond_distance ** 2
    sites = structure.sites
    num_sites = len(sites)
    if num_sites < 2:
        return []
    
    prop = []
    for site in sites:
        prop.append({
            'element': site.specie.symbol,
            'electroneg': site.specie.X,
            'is_metal': site.specie.is_metal,
            'is_nonmetal': site.specie.is_nonmetal,
            'covalent_radius': covalent_radii[site.specie.Z]
        })
    
    max_radius = max(covalent_radii)
    max_cutoff = max_radius * 2 * max_distance_ratio # Maximum distance to consider bonding
    spatial, cell_size = setup_spatial_grid(sites, max_cutoff)
    bonds = []
    closet = {}

    for i in range(num_sites):
        xi, yi, zi = sites[i].coords
        pi = prop[i]

        for j in get_neighboring_candidates(sites[i].xyz, sites, spatial, cell_size):
            if j <= i:
                continue
            xj, yj, zj = sites[j].xyz
            pj = prop[j]
            dx, dy, dz = xi - xj, yi - yj, zi - zj
            dist_sq = dx * dx + dy * dy + dz * dz
            dist = m.sqrt(dist_sq)
            if dist_sq < min_dist_sq or not pi["covalent_radius"] or not pj["covalent_radius"]:
                continue
            expected = pi["covalent_radius"] + pj["covalent_radius"]

            if dist > expected * max_distance_ratio:
                continue

            en_diff = abs(pi["electroneg"] - pj["electroneg"])
            en_ratio = en_diff / (pi["electroneg"] + pj["electroneg"])

            bond_strength = 1.
            if pi["is_metal"] and pj["is_metal"]:
                bond_strength *= metal_metal_penalty
            elif (pi["is_metal"] and pj["is_nonmetal"]) or (pi["is_nonmetal"] and pj["is_metal"]):
                bond_strength *= metal_nonmetal_bonus
            if en_diff > electroneg_threshold:
                bond_strength *= 1.3
            elif en_diff < 0.5:
                bond_strength *= similarity_x_bonus

            dist_weight = m.exp(-(dist/expected - 1) **2 / 0.18)
            en_weight = 1.0 - 0.3 * en_ratio
            strength = bond_strength * dist_weight * en_weight

            if pi["element"] == pj["element"]:
                strength *= same_element_penalty

            ca = closet.get(i, float('inf'))
            cb = closet.get(j, float('inf'))

            if dist > ca:
                strength *= m.exp(-(dist/ca - 1)/0.5)
            if dist > cb:
                strength *= m.exp(-(dist/cb - 1)/0.5)
            
            if strength > strength_threshold:
                bonds.append({"pos_1": sites[i].xyz,
                              "pos_2": sites[j].xyz,
                              "index_1": i,
                              "index_2": j,
                              "bond_length": dist,
                              "strength": strength,
                              "transform_matrix": compute_bond_transform(sites[i].xyz, sites[j].xyz)})
    return bonds
    

def main():
    Structure.from_file("working/surface_model/surface/miller_010_layers2_minlat10.0_symTrue/slab_1.cif")

if __name__ == '__main__':
    main()
    




