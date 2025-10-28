#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

import argparse
import random

from ase import Atoms
from ase.io import read, write
from pymatgen.core import Species

from build_model import generate_water_mol_box_packmol


def build_water_box(xyz_file: str, ions: list, ions_num: list):
    """
    Based on the input xyz file, build a water box with specified ions and their numbers.
    """
    ...


def build_water_box(a, b, c, ions, ions_num):
    """
    Build a water box with specified dimensions and ions.
    """
    ...


def parser_args():

    parser = argparse.ArgumentParser(
        description="Build a water box with specified ions"
    )
    parser.add_argument("-x", "--xyz_file", type=str, help="Path to the input xyz file")
    parser.add_argument("-a", type=float, help="Box length in x direction")
    parser.add_argument("-b", type=float, help="Box length in y direction")
    parser.add_argument("-c", type=float, help="Box length in z direction")
    parser.add_argument(
        "-s", "--num_sol", type=int, default=192, help="Number of water molecules"
    )
    parser.add_argument(
        "-i", "--ions", nargs="+", default=None, help="List of ion types to include"
    )
    parser.add_argument(
        "-n",
        "--ions_num",
        nargs="+",
        default=None,
        type=int,
        help="List of ion counts corresponding to ion types",
    )
    args = parser.parse_args()
    return args


def main():
    args = parser_args()

    ions = args.ions if args.ions else ["K+", "Cl-"]
    ions_number = args.ions_num if args.ions_num else [1, 1]
    num_sol = args.num_sol

    if ions and ions_number:
        if len(ions) != len(ions_number):
            raise ValueError("The length of ions and ions_number must be the same.")

    charge = 0
    extra_mol = []
    extra_num = []
    for ion, ion_num in zip(ions, ions_number):
        ion_s = Species(ion)
        if ion_s.oxi_state:
            charge += ion_s.oxi_state * ion_num
        else:
            raise ValueError(f"The ion {ion} does not have a valid oxidation state.")
        extra_mol.append(Atoms(str(ion_s.element), positions=[[0, 0, 0]]))
        extra_num.append(ion_num)
        # print(extra_mol)

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
        print("Warning: ions and ions_number must be both provided.")
    else:
        msg = "No extra ions added."

    a, b, c = (
        read(args.xyz_file, index=0).get_cell_lengths_and_angles()[:3]
        if args.xyz_file
        else (args.a, args.b, args.c)
    )
    solvation = "H2O"
    water_box = generate_water_mol_box_packmol(
        a,
        b,
        c,
        sol=solvation,
        num=num_sol,
        molecule=extra_mol,
        molecule_number=extra_num,
        verbose=False,
        region="top",
        seed=random.randint(1, 10000),
    )
    write(
        f"water_box_{[ion for ion in ions]}_{[num for num in ions_number]}.xyz",
        water_box,
    )


if __name__ == "__main__":
    main()
