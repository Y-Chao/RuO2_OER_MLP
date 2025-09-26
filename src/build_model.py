#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"


"""
Build interface RuO2 models, including surface slab and solvation/ions.
Functions:
- load_bulk_db: Load bulk structures from ASE database.
- build_surface_model: Build surface slab from bulk structure.
- build_interface_model: Add solvation and ions to the surface slab.
"""

import argparse
import os

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.db import connect
from ase.io import read, write
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor


def load_bulk_db(db_name: str):
    """
    Load the database, which store the bulk infomation with the following keys: [sample, crystal, xc]
    Among them,
        `sample` represents the composition of bulk structure,
        `crystal` represents the crystal shape, like rutile, fcc
        `xc` represnets the functional used for bulk oprimization
    """
    db = connect(db_name)
    return {
        row.sample + "-" + row.crystal + "-" + row.xc: db.get_atoms(row.id)
        for row in db.select()
    }


def check_termination(slab: Slab, terminations: list[str] | str) -> bool:
    """
    Check if the slab has the desired terminations.
    terminations: list of str, e.g., ["O", "Ru"]
    """
    top_layer = slab.cart_coords[:, 2].max()
    top_index = np.where(slab.cart_coords[:, 2] > top_layer - 0.3)[
        0
    ]  # 0.3 is a threshold for layer

    top_elements = set(slab.species[i].symbol for i in top_index)
    print(f"Top layer elements: {top_elements}")

    if isinstance(terminations, str):
        terminations = [terminations]
    return set(terminations) == top_elements


def move_to_bottom(slab: Slab | Atoms) -> Atoms:
    """
    Move the slab to the bottom of the cell.
    """
    if isinstance(slab, Slab):
        slab = AseAtomsAdaptor.get_atoms(slab)
    z_min = slab.positions[:, 2].min()
    slab.positions[:, 2] -= z_min
    slab.wrap()
    return slab


def extend_vacuum(slab: Slab, vacuum: float, orig_vacuum: float) -> Slab:
    """
    Extend the vacuum of the slab to the desired value.
    Structure in pymatgen is hard to modify, so we first convert it to ASE Atoms,
    then modify the lattice.
    """
    ase_slab = AseAtomsAdaptor.get_atoms(slab)
    new_lattice = np.array(ase_slab.cell)
    new_lattice[2, 2] += vacuum - orig_vacuum
    ase_slab.set_cell(new_lattice)
    ase_slab = move_to_bottom(ase_slab)
    return ase_slab


def make_simple_supercell(slab: Atoms | Slab, min_lattice: float) -> Atoms:
    """
    Make the slab into a supercell with minimum lattice constant.
    """
    if isinstance(slab, Slab):
        slab = AseAtomsAdaptor.get_atoms(slab)
    a, b = slab.cell.cellpar()[:2]
    na = int(np.ceil(min_lattice / a))
    nb = int(np.ceil(min_lattice / b))
    return slab.repeat((na, nb, 1))


def supercell_transformation(
    slab: Atoms | Slab, transition_matrix: np.ndarray
) -> Atoms:
    """
    This function is referred from pymatgen's SupercellTransformation class.

    Make the slab into a supercell with transition matrix.
    transition_matrix: 2D array, e.g., [[1, 0], [0, 1]] means no change,
    [[2, 0], [0, 1]] means double the a lattice, [[1, 1], [0, 1]] means
    make a + b as new a lattice.
    """
    if isinstance(slab, Slab):
        slab = AseAtomsAdaptor.get_atoms(slab)

    new_a, new_b = redefine_lattice(a, b, transition_matrix)
    na = int(np.round(new_a / a))
    nb = int(np.round(new_b / b))
    return slab.repeat((na, nb, 1))


def make_supercell(slab: Atoms, transition_matrix: np.ndarray) -> Atoms:
    """
    Make the slab into a supercell with transition matrix.
    transition_matrix: 2D array, e.g., [[1, 0], [0, 1]] means no change,
    [[2, 0], [0, 1]] means double the a lattice, [[1, 1], [0, 1]] means
    make a + b as new a lattice.
    """
    is_int = np.isclose(transition_matrix, np.round(transition_matrix))
    if not is_int.all():
        raise ValueError("Transition matrix must be integer.")

    # Ensure the TM is a positive definite matrix
    if not np.all(np.linalg.eigvals(transition_matrix) > 0):
        raise ValueError("Transition matrix must be positive definite.")

    def redefine_lattice(a, b, tm):
        new_a = tm[0, 0] * a + tm[0, 1] * b
        new_b = tm[1, 0] * a + tm[1, 1] * b
        return new_a, new_b


def calc_sol_num(
    slab: Atoms, sol: Atoms, c_sup=None, density=1.0, surf_height=2.0, scale=1.0
):
    """Calculate the number of solvent molecules according to the density of the solvent.
    Parameters:
    -----------
    slab: ase.Atoms
        The slab surface
    sol: ase.Atoms
        The solution molecular
    c_sup: float
        The c axis of the supercell, default is None.
    density: float
        The density of the solvent, default is 1.0 g/mL
    surf_height: float
        The height of above the slab, default is 2.0 Ang
    Returns:
    --------
    num: int
        The number of solvent molecules
    """

    (a, b, c) = slab.cell.array
    assert np.allclose(
        slab.cell.cellpar()[3:5], 90.0
    ), "Currently, it only supports orthogonal slabs."

    slab_area = np.linalg.norm(np.cross(a, b))  # units: Angstrom^2
    if c_sup is not None:
        z = c_sup  # units: Angstrom
    else:
        # Get the maximum z height of the metal atoms
        metal_z_height = [atom.position[2] for atom in slab if atom.number > 18]
        z = np.linalg.norm(c) - max(metal_z_height) - surf_height  # units: Angstrom
    sol_volume = slab_area * z * 10**-24  # units: mL
    sol_mass = np.sum(sol.get_masses())  # units: g/mol
    num = int(sol_volume * scale * density / sol_mass * mol)
    print("The number of solvent molecules is: {}".format(num))
    return num


def generate_water_box_packmol(a, b, c, num, output_file="water_box.pdb"):
    """
    Generate a water box using packmol.
    Parameters:
    -----------
    a, b, c: float
        The lattice constants of the box.
    num: int
        The number of water molecules.
    output_file: str
        The output file name.
    Returns:
    --------
    water_box: ase.Atoms
        The water box.
    """
    input_file = "water_box.in"
    with open(input_file, "w") as f:
        f.write("# The input file for solvation box build by packmol\n")
        f.write("tolerance 2.0\n")
        f.write("output {}\n".format(output_file))
        f.write("structure sol.pdb\n")
        f.write("  number {}\n".format(num))
        f.write("  inside box 0. 0. 0.  {} {} {}\n".format(a - 1.0, b - 1.0, c - 2.0))
        f.write("end structure\n")

    exit_status = os.system("packmol < {}".format(input_file))
    if exit_status != 0:
        print("Error in conducting packmol. Please check the packmol installation.")
        for f in ["sol.pdb", "water_box.pdb", "water_box.in"]:
            if os.path.exists(f):
                os.remove(f)
        return None
    water_box = read("water_box.pdb")
    for f in ["sol.pdb", "water_box.pdb", "water_box.in"]:
        if os.path.exists(f):
            os.remove(f)
    return water_box


def generate_water_mol_box_packmol(a, b, c, num, output_file="water_mol.pdb", **kwargs):
    """
    Generate a water box with other molecules using packmol.
    Parameters:
    -----------
    a, b, c: float
        The lattice constants of the box.
    num: int
        The number of water molecules.
    output_file: str
        The output file name.
    molecule: ase.Atoms
        The other molecule to be added.
    molecule_number: int
        The number of other molecules.
    molecule_file: str
        The file name of the other molecule.
    Returns:
    --------
    water_box: ase.Atoms
        The water box with other molecules.
    """
    input_file = "water_box.in"
    mol = kwargs.get("molecule", None)
    mol_num = kwargs.get("molecule_number", None)
    mol_file = kwargs.get("molecule_file", None)

    with open(input_file, "w") as f:
        f.write("# The input file for solvation box build by packmol\n")
        f.write("tolerance 2.0\n")
        f.write("output {}\n".format(output_file))
        f.write("structure sol.pdb\n")
        f.write("  number {}\n".format(num))
        f.write("  inside box 0. 0. 0.  {} {} {}\n".format(a - 1.0, b - 1.0, c - 2.0))
        f.write("end structure\n")
        if mol_file is not None:
            f.write("structure mol.pdb\n")
            f.write("  number {}\n".format(mol_num))
            f.write(
                "  inside box 0. 0. 0.  {} {} {}\n".format(a - 1.0, b - 1.0, c - 2.0)
            )
            f.write("end structure\n")
        else:
            if mol is not None:
                write("mol.pdb", mol)
                f.write("structure mol.pdb\n")
                f.write("  number {}\n".format(mol_num))
                f.write(
                    "  inside box 0. 0. 0.  {} {} {}\n".format(
                        a - 1.0, b - 1.0, c - 2.0
                    )
                )
                f.write("end structure\n")

    exit_status = os.system("packmol < {}".format(input_file))

    if exit_status != 0:
        print("Error in conducting packmol. Please check the packmol installation.")
        for f in ["sol.pdb", "water_box.pdb", "water_box.in"]:
            if os.path.exists(f):
                os.remove(f)
        return None
    water_box = read("water_box.pdb")
    for f in ["sol.pdb", "water_box.pdb", "water_box.in"]:
        if os.path.exists(f):
            os.remove(f)
    return water_box


def add_water_box(slab, water_box, surf_height):
    """Add the water box above the slab surface.
    Parameters:
    -----------
    slab: ase.Atoms
        The slab surface
    water_box: ase.Atoms
        The water box
    Returns:
    --------
    slab_water: ase.Atoms
        The slab with water box
    """
    z_max = slab.positions[:, 2].max()
    water_box.positions[:, 2] += z_max + surf_height
    slab_water = slab + water_box
    slab_water.wrap()
    return slab_water


def build_surface_model(
    bulk, miller_index, layers, vacuum, terminations=None, symetry=True
):
    bulk_pmg = AseAtomsAdaptor.get_structure(bulk)
    slab_gen = SlabGenerator(
        initial_structure=bulk_pmg,
        miller_index=miller_index,
        min_slab_size=layers,
        min_vacuum_size=1,  # Fix vacuum to 1 layer, and later adjust it based on the lattice c
        lll_reduce=False,  # lll_reduce means to make the slab orthogonal, which is not desired here
        in_unit_planes=True,
        primitive=True,
        max_normal_search=20,
        reorient_lattice=True,
    )
    vacuum_size = slab_gen.oriented_unit_cell.lattice.c
    if vacuum_size > 10:  # Avoid too large OUC
        raise ValueError(
            f"The vacuum size is too large: {vacuum_size}. Please check the miller index and layers."
        )
    tmp_surfaces = slab_gen.get_slabs(
        ftol=0.5, symmetrize=symetry, filter_out_sym_slabs=True
    )

    if len(tmp_surfaces) == 0:
        raise ValueError("No surface generated. Please check the input parameters.")

    # Check the terminations
    surfaces = []
    if terminations:
        for s in tmp_surfaces:
            if check_termination(s, terminations):
                surfaces.append(s)
    else:
        surfaces = tmp_surfaces

    # Apply the define vacuum size
    for i in range(len(surfaces)):
        surfaces[i] = extend_vacuum(surfaces[i], vacuum, vacuum_size)
    return surfaces


def build_interface_model(
    surface: Atoms | Slab,
    solvation="H2O",
    pH: str = "neutral",
    ions: list[str] | None = None,
    ions_number: list[int] | None = None,
):
    """
    Add solvation and ions to the surface slab.
    solvation: str, e.g., "H2O"
    pH: str: "acidic", "neutral", "basic"
    ions: list of str, e.g., ["Na", "Cl"]

    To be implemented:
    1. pH effect
        a. acidic: add H3O+ ions
        b. basic: add OH- ions and K+ ions
        c. neutral: add H2O
    2. ionic number
    """
    if isinstance(surface, Slab):
        surface = AseAtomsAdaptor.get_atoms(surface)

    write("sol.pdb", molecule("H2O"))
    num_sol = calc_sol_num(
        surface, molecule("H2O"), density=1.0, surf_height=2.0, scale=1.0
    )

    # Handle pH and ions
    extra_mol = None
    extra_num = None
    extra_file = None

    if pH == "acidic":
        try:
            extra_mol = molecule("H3O")
        except Exception:
            extra_mol = molecule("H3O+")
        extra_num = 1
    elif pH == "basic":
        try:
            extra_mol = molecule("OH")
        except Exception:
            extra_mol = molecule("OH-")
        extra_num = 1
    elif ions:
        try:
            extra_mol = molecule(ions[0])
            extra_num = 1
        except Exception:
            extra_mol = None
            extra_num = None

    if extra_mol is not None and extra_num is not None:
        water_box = generate_water_mol_box_packmol(
            *surface.cell.cellpar()[:3],
            num=num_sol,
            molecule=extra_mol,
            molecule_number=extra_num,
        )
        msg = f"Interface model with {num_sol} water molecules and extra species generated."
    else:
        water_box = generate_water_box_packmol(*surface.cell.cellpar()[:3], num=num_sol)
        msg = f"Interface model with {num_sol} water molecules generated."

    if water_box is None:
        raise RuntimeError("Failed to generate water box.")

    interface = add_water_box(surface, water_box, surf_height=2.0)
    write("interface.xyz", interface)
    print(msg)
    return interface


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Build RuO2 interface models",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-b",
        "--bulk_db",
        type=str,
        default="RuO2_bulk.db",
        help="Path to the bulk database",
    )
    parser.add_argument(
        "--crystal",
        type=str,
        default="RuO2-rutile-pbe",
        help="Crystal structure of RuO2",
    )
    parser.add_argument(
        "--miller",
        nargs=3,
        type=int,
        default=[1, 0, 0],
        help="Miller index for surface generation",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    bulks_results = load_bulk_db(args.bulk_db)
    print(f"Loaded {len(bulks_results)} bulk structures from {args.bulk_db}")
    print(f"Available structures: {list(bulks_results.keys())}")


if __name__ == "__main__":
    main()
