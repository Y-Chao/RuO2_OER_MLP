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

To be implemented:
- add water molecules on non-orthogonal slabs.
"""

import argparse
import subprocess

import numpy as np
from ase import Atoms, build
from ase.constraints import FixAtoms
from ase.db import connect
from ase.io import read, write
from ase.units import mol
from monty.string import boxed, indent, marquee
from monty.tempfile import ScratchDir
from pymatgen.core.periodic_table import Species
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

from load_params import load_surface_toml


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


def fix_surface(slab: Atoms | Slab, fix: str, layers: int) -> Atoms:
    """
    Fix the slab.
    Parameters:
    -----------
    slab: ase.Atoms or pymatgen Slab object
    fix: str, "top", "bottom", "all", "none"
    layers: int, number of layers to fix
    Returns:
    --------
    slab: ase.Atoms
    """
    if isinstance(slab, Slab):
        slab = AseAtomsAdaptor.get_atoms(slab)

    if fix == "none":
        return slab
    elif fix == "bottom":
        fix_layer = np.floor(layers / 2)
        slab_height = slab.positions[:, 2].max() - slab.positions[:, 2].min()
        layer_height = slab_height / layers
        z_cutoff = slab.positions[:, 2].min() + fix_layer * layer_height
        fix_indices = []
        for atom in slab:
            if atom.position[2] < z_cutoff:
                fix_indices.append(atom.index)
        slab.set_constraint(FixAtoms(indices=fix_indices))
        return slab

    elif fix == "center":
        z_center = (slab.positions[:, 2].max() + slab.positions[:, 2].min()) / 2
        fix_layer = np.floor(layers / 2)
        slab_height = slab.positions[:, 2].max() - slab.positions[:, 2].min()
        layer_height = slab_height / layers
        z_lower = z_center - (fix_layer / 2) * layer_height
        z_upper = z_center + (fix_layer / 2) * layer_height
        fix_indices = []
        for atom in slab:
            if atom.position[2] < z_lower or atom.position[2] > z_upper:
                fix_indices.append(atom.index)
        slab.set_constraint(FixAtoms(indices=fix_indices))
        return slab
    else:
        raise ValueError(
            f"Invalid fix option: {fix}. Choose from 'top', 'bottom', 'all', 'none'."
        )


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
    # print(f"Top layer elements: {top_elements}")

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
    slab: Atoms, sol: Atoms | str, c_sup=None, density=1.0, surf_height=2.0, scale=1.0
):
    """Calculate the number of solvent molecules according to the density of the solvent.
    Parameters:
    -----------
    slab: ase.Atoms
        The slab surface
    sol: ase.Atoms | str
        The solution molecular or the name of the solvent molecular.
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

    # Slab section
    (a, b, c) = slab.cell.array
    assert np.allclose(
        slab.cell.cellpar()[3:5], 90.0, rtol=0.5
    ), "Currently, it only supports orthogonal slabs."

    slab_area = np.linalg.norm(np.cross(a, b))  # units: Angstrom^2
    if c_sup is not None:
        z = c_sup  # units: Angstrom
    else:
        # Get the maximum z height of the metal atoms
        metal_z_height = [atom.position[2] for atom in slab if atom.number > 18]
        z = np.linalg.norm(c) - max(metal_z_height) - surf_height  # units: Angstrom
    sol_volume = slab_area * z * 10**-24  # units: mL

    # Solvation section
    if isinstance(sol, str):
        sol = build.molecule(sol)
    sol_mass = np.sum(sol.get_masses())  # units: g/mol
    num = int(sol_volume * scale * density / sol_mass * mol)
    print("The number of solvent molecules is: {}".format(num))
    return num


def generate_water_box_packmol(
    a, b, c, sol, num, output_file="water_box.pdb", verbose=False, **kwargs
):
    """
    Generate a water box using packmol.
    Parameters:
    -----------
    a, b, c: float
        The lattice constants of the box.
    sol: ase.Atoms | str
        The solvent molecular or the name of the solvent molecular.
    num: int
        The number of water molecules.
    output_file: str
        The output file name.
    verbose: bool
        Whether to keep the temporary files.
    z_height: float
        The height of above the slab, default is 1.0 Angstrom.
        This will keep the water molecules away from the both sides of slab surface 1.0 Angstrom.
    Returns:
    --------
    water_box: ase.Atoms
        The water box.
    """
    input_file = "water_box.in"
    if verbose:
        copy_to_current = True
    else:
        copy_to_current = False

    z_height = kwargs.get("surface_height", 1.0)

    with ScratchDir(
        rootpath=".",
        copy_to_current_on_exit=copy_to_current,
        delete_removed_files=False,
    ) as temp:
        if isinstance(sol, str):
            sol = build.molecule(sol)
        write("sol.pdb", sol)
        with open(input_file, "w") as f:
            f.write("# The input file for solvation box build by packmol\n")
            f.write("tolerance 2.0\n")
            f.write("output {}\n".format(output_file))
            f.write("structure sol.pdb\n")
            f.write("  number {}\n".format(num))
            f.write(
                "  inside box 0. 0. 0.  {} {} {}\n".format(
                    a - 1.0, b - 1.0, c - z_height * 2.0
                )
            )
            f.write("end structure\n")

        with open(input_file, "r") as f:
            if verbose:
                out = None
            else:
                out = subprocess.DEVNULL
            results = subprocess.run(
                ["packmol"],
                stdin=f,
                stdout=out,
                # stderr=out,
            )

        # exit_status = os.system("packmol < {}".format(input_file))
        if results.returncode != 0:
            print("Error in conducting packmol. Please check the packmol installation.")

        water_box = read("water_mol.pdb")
    return water_box


def generate_water_mol_box_packmol(
    a, b, c, sol, num, output_file="water_mol.pdb", verbose=False, **kwargs
):
    """
    Generate a water box with other molecules using packmol.
    Parameters:
    -----------
    a, b, c: float
        The lattice constants of the box.
    sol: ase.Atoms | str
        The solvent molecular or the name of the solvent molecular.
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
    verbose: bool
        Whether to keep the temporary files.
    z_height: float
        The height of above the slab, default is 1.0 Angstrom.
        This will keep the water molecules away from the both sides of slab surface 1.0 Angstrom.
    region: str
        The region to place the other molecules, default is "bottom".
    Returns:
    --------
    water_box: ase.Atoms
        The water box with other molecules.
    """
    input_file = "water_box.in"
    mol = kwargs.get("molecule", None)
    mol_num = kwargs.get("molecule_number", None)
    region = kwargs.get("region", "bottom")  # Default to bottom region
    z_region = {
        "bottom": c / (2 * 3),
        "middle": c / 3,
        "top": c / 2,
    }  # [bottom, lower, middle, top]
    z_height = kwargs.get("surface_height", 1.0)

    if verbose:
        copy_to_current = True
    else:
        copy_to_current = False

    with ScratchDir(
        rootpath=".",
        copy_to_current_on_exit=copy_to_current,
        delete_removed_files=False,
    ) as temp:
        if isinstance(sol, str):
            sol = build.molecule(sol)
        write("sol.pdb", sol)
        with open(input_file, "w") as f:
            f.write("# The input file for solvation box build by packmol\n")
            f.write("tolerance 2.0\n")
            f.write("output {}\n".format(output_file))
            f.write("structure sol.pdb\n")
            f.write("  number {}\n".format(num))
            f.write(
                "  inside box 0. 0. 0.  {} {} {}\n".format(
                    a - 1.0, b - 1.0, c - z_height * 2.0
                )
            )
            f.write("end structure\n")
            if mol is not None and len(mol) == len(mol_num):
                for i in range(len(mol)):
                    write(f"mol_{i}.pdb", mol[i])
                    f.write(f"structure mol_{i}.pdb\n")
                    f.write("  number {}\n".format(mol_num[i]))
                    f.write(
                        "  inside box 0. 0. 0.  {} {} {}\n".format(
                            a - 1.0, b - 1.0, z_region[region]
                        )
                    )
                    f.write("end structure\n")

        with open(input_file, "r") as f:
            if verbose:
                out = None
            else:
                out = subprocess.DEVNULL
            results = subprocess.run(
                ["packmol"],
                stdin=f,
                stdout=out,
                # stderr=out,
            )
        # exit_status = os.system("packmol < {}".format(input_file))

        if results.returncode != 0:
            print("Error in conducting packmol. Please check the packmol installation.")

        water_box = read("water_mol.pdb")
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
    bulk,
    miller_index,
    layers,
    vacuum,
    symmetry,
    terminations=None,
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
        ftol=0.5, symmetrize=symmetry, filter_out_sym_slabs=True
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
    region: str = "bottom",
    verbose: bool = False,
    **kwargs,
):
    """
    Add solvation and ions to the surface slab.
    solvation: str, e.g., "H2O"
    pH: str: "acidic", "neutral", "basic"
    ions: list of str, e.g., ["Na+", "Cl-"]

    To be implemented:
    1. pH effect
        a. acidic: add H3O+ ions
        b. basic: add OH- ions and K+ ions
        c. neutral: add H2O
    2. ionic number
    """
    if isinstance(surface, Slab):
        surface = AseAtomsAdaptor.get_atoms(surface)

    num_sol = calc_sol_num(
        surface, build.molecule("H2O"), density=1.0, surf_height=2.0, scale=1.0
    )

    # Handle pH and ions
    extra_mol = []
    extra_num = []
    # extra_file = None

    # First handle the ions effect, then pH effect
    # based on the ions number, and equilibrium
    if ions and ions_number:
        if len(ions) != len(ions_number):
            raise ValueError("The length of ions and ions_number must be the same.")

        charge = 0
        for ion, ion_num in zip(ions, ions_number):
            ion_s = Species(ion)
            if ion_s.oxi_state:
                charge += ion_s.oxi_state * ion_num
            else:
                raise ValueError(
                    f"The ion {ion} does not have a valid oxidation state."
                )
            extra_mol.append(Atoms(str(ion_s.element), positions=[[0, 0, 0]]))
            extra_num.append(ion_num)

        if charge > 0 and pH == "acidic":
            msg = f"Adding H3O+ is contradictory to the positive ions {ions}."
            raise ValueError(msg)
        elif charge < 0 and pH == "basic":
            msg = f"Adding OH- is contradictory to the negative ions {ions}."
            raise ValueError(msg)
        elif charge != 0 and pH == "neutral":
            msg = "Charge imbalance without H3O+."
            raise ValueError(msg)

        if charge > 0:
            extra_mol += [Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.9742]])]
            extra_num += [int(charge)]

        elif charge < 0:
            extra_mol += [
                Atoms(
                    "H3O",
                    positions=[
                        [0.0, 0.9441, -0.2125],
                        [0.8176, -0.4720, -0.2125],
                        [-0.8176, -0.4720, -0.2125],
                        [0.0, 0.0, 0.0797],
                    ],
                )
            ]
            extra_num += [int(-charge)]
    elif ions or ions_number:
        raise ValueError("Both ions and ions_number must be provided.")
    else:
        msg = "No extra ions added."

    slab_cell = surface.cell.cellpar()
    a, b, c = slab_cell[:3]
    c_sup = c - surface.positions[:, 2].max()

    if extra_mol is not None and extra_num is not None:
        water_box = generate_water_mol_box_packmol(
            a,
            b,
            c_sup,
            sol=solvation,
            num=num_sol,
            molecule=extra_mol,
            molecule_number=extra_num,
            verbose=verbose,
            region=region,
        )
        msg = f"Interface model with {num_sol} water molecules and extra species {','.join([str(e_mol.get_chemical_formula()) for e_mol in extra_mol])} generated."
    else:
        water_box = generate_water_box_packmol(a, b, c_sup, sol=solvation, num=num_sol)
        msg = f"Interface model with {num_sol} water molecules generated."

    if water_box is None:
        raise RuntimeError("Failed to generate water box.")

    interface = add_water_box(
        surface, water_box, surf_height=kwargs.get("surface_height", 1.0)
    )
    # write("interface.xyz", interface)
    print(msg)
    return interface


def detailed_bulk(bulk: Atoms, args) -> str:
    """
    Get the detailed information of the bulk structure.
    """
    formula = bulk.get_chemical_formula()
    num_atoms = len(bulk)
    cell = bulk.cell.cellpar()
    a, b, c, alpha, beta, gamma = cell
    info = [
        f"Bulk formula: ({formula})",
        f"Num atoms: {num_atoms}",
        f"Crystal: {args.bulk_info.crystal}",
        f"Functional: {args.bulk_info.xc}",
        f"abc   : {a:.2f}, {b:.2f}, {c:.2f}",
        f"alpha beta gamma: {alpha:.2f}, {beta:.2f}, {gamma:.2f}",
    ]
    return "\n".join(info)


def detailed_slab(slab: Atoms, args) -> str:
    """
    Get the detailed information of the slab.
    """
    formula = slab.get_chemical_formula()
    num_atoms = len(slab)
    cell = slab.cell.cellpar()
    a, b, c, alpha, beta, gamma = cell
    info = [
        f"Slab formula: ({formula})",
        f"Num atoms: {num_atoms}",
        f"Miller index: ({args.slab_info.miller_index})",
        f"abc   : {a:.2f}, {b:.2f}, {c:.2f}",
        f"alpha beta gamma: {alpha:.2f}, {beta:.2f}, {gamma:.2f}",
        f"Vacuum: {args.slab_info.vacuum} Angstrom",
        f"Layers: {args.slab_info.layers}",
        f"Terminations: {args.slab_info.terminations}",
        f"Symmetry: {args.slab_info.symmetry}",
        f"Fix: {args.slab_info.fix}",
        f"Min lattice: {args.slab_info.min_lattice} Angstrom",
    ]
    return "\n".join(info)


def detailed_interface(interface: Atoms, args) -> str:
    """
    Get the detailed information of the interface.
    """
    formula = interface.get_chemical_formula()
    num_atoms = len(interface)
    cell = interface.cell.cellpar()
    a, b, c, alpha, beta, gamma = cell
    info = [
        f"Interface formula: ({formula})",
        f"Num atoms: {num_atoms}",
        f"abc   : {a:.2f}, {b:.2f}, {c:.2f}",
        f"alpha beta gamma: {alpha:.2f}, {beta:.2f}, {gamma:.2f}",
        f"Solvation: {args.interface_info.solvation}",
        f"pH: {args.interface_info.pH}",
        f"Ions: {args.interface_info.ions}",
        f"Ions number: {args.interface_info.ions_number}",
        f"Surface height: {args.interface_info.surface_height} Angstrom",
        f"Verbose: {args.interface_info.verbose}",
    ]
    return "\n".join(info)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Build RuO2 interface models",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--tomp_file",
        type=str,
        default="xMOI_configurations.toml",
        help="Path to the TOML configuration file",
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


def run():
    args = parse_arguments()
    package_info = boxed("Extended Metal Oxide Interface", pad=5)
    print(package_info)
    print("\n")

    ############################ Load configuratiuons ############################
    conf = load_surface_toml(args.tomp_file)

    ############################ Load bulk structures ############################
    task_info = marquee("Load bulk structures", width=78)
    print(task_info)
    bulks_results = load_bulk_db(conf.bulk_info.bulk_db)
    print(f"Loaded {len(bulks_results)} bulk structures from {conf.bulk_info.bulk_db}")
    intent_info = indent(f"Available structures: {list(bulks_results.keys())}", 4)
    print(intent_info)
    # pprint(list(bulks_results.keys()))
    print("\n")
    print(indent(f"Selected sample: {conf.bulk_info.sample}", 4))
    print(indent(f"Selected bulk structure: {conf.bulk_info.crystal}", 4))
    print(indent(f"Functional: {conf.bulk_info.xc}", 4))
    print("\n")

    sample_name = (
        conf.bulk_info.sample + "-" + conf.bulk_info.crystal + "-" + conf.bulk_info.xc
    )
    if sample_name not in bulks_results:
        raise ValueError(
            f"The bulk structure {sample_name} is not found in the database."
        )

    bulk = bulks_results[sample_name]
    print("Bulk structure details:")
    print(indent(f"{detailed_bulk(bulk, conf)}", 4))
    print("\n")

    ########################## Generate surface models ##########################
    task_info = marquee("Generate surface models", width=78)
    print(task_info)
    print(indent(f"Miller index: {conf.slab_info.miller_index}", 4))
    print(indent(f"Layers: {conf.slab_info.layers}", 4))
    print(indent(f"Vacuum: {conf.slab_info.vacuum}", 4))
    print(indent(f"Terminations: {conf.slab_info.terminations}", 4))
    print(indent(f"Symmetry: {conf.slab_info.symmetry}", 4))

    surfaces = build_surface_model(
        bulk,
        miller_index=conf.slab_info.miller_index,
        layers=conf.slab_info.layers,
        vacuum=conf.slab_info.vacuum,
        terminations=conf.slab_info.terminations,
        symmetry=conf.slab_info.symmetry,
    )

    print(f"Generated {len(surfaces)} surface models.")
    for i, s in enumerate(surfaces):
        slab = fix_surface(
            s,
            fix=conf.slab_info.fix,
            layers=conf.slab_info.layers,
        )
        slab = make_simple_supercell(slab, conf.slab_info.min_lattice)
        # write(f"slab_{i+1}.xyz", slab)
        print(f"Surface model {i+1}")
        print(indent(f"{detailed_slab(slab, conf)}", 4))
        if conf.interface_info.solvation is None:
            write(f"slab_{i+1}.xyz", slab)
        print("\n")

    ########################## Generate interface models ##########################
    if conf.interface_info.solvation:
        task_info = marquee("Generate interface models", width=78)
        print(task_info)
        print(indent(f"Solvation: {conf.interface_info.solvation}", 4))
        print(indent(f"pH: {conf.interface_info.pH}", 4))
        print(indent(f"Ions: {conf.interface_info.ions}", 4))
        print(indent(f"Ions number: {conf.interface_info.ions_number}", 4))
        print("\n")
        for i, s in enumerate(surfaces):
            interface = build_interface_model(
                slab,
                solvation=conf.interface_info.solvation,
                surface_height=conf.interface_info.surface_height,
                pH=conf.interface_info.pH,
                ions=conf.interface_info.ions,
                ions_number=conf.interface_info.ions_number,
                region="bottom",
                verbose=conf.interface_info.verbose,
            )
            write(f"interface_{i+1}.xyz", interface)
            print(f"Interface model {i+1}")
            print(indent(f"{detailed_interface(interface, conf)}", 4))
            print("\n")


if __name__ == "__main__":
    run()
