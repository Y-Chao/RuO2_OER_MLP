#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = "Chao Yang"
__version__ = "1.0"

import argparse
import glob
import os

import numpy as np
from ase.calculators.vasp import Vasp
from ase.db import connect
from ase.io import read
from monty.string import indent, marquee
from mp_api.client import MPRester

"""
Build the bulk database from materials project.
Targets:
    - Download the needed bulk structures from materials project.
    - Store them in a local ase.db database.
    - Prepare the calculation input files.
    - Perform the DFT calculations.
"""


def query_structure(chemsystem, fields=None):
    """
    Query the bulk structures from materials project.
    Parameters
    ----------
    chemsystem : list
        List of chemical systems to query.
    fields : list, optional
        List of fields to query, by default None.
    Returns
    -------
    docs : list
        List of queried documents.
    """
    api_mpg = os.environ.get("mp_api", None)
    if api_mpg is None:
        raise ValueError("Please set the environment variable mp_api")

    with MPRester(api_mpg) as mpr:
        docs = mpr.get_entries(
            chemsys_formula_mpids=chemsystem,
            property_data=fields,
            conventional_unit_cell=True,
        )
    return docs


def parse_bulk_data(dbname, docs, xc="GGA"):
    """
    Store the queried bulk structures in a local ase.db database.
    Parameters
    ----------
    dbname : str
        Name of the database.
    docs : list
        List of queried documents.
    xc : str, optional
        Exchange-correlation functional, by default "GGA".
    Returns
    -------
    db : ase.db.connection
        Connection to the local database.
    """
    with connect(dbname) as db:
        MPID = []
        for doc in docs:
            # Filter the documents based on the xc and avoid duplicates
            if doc.data["run_type"] != xc:
                continue
            if doc.data["material_id"] in MPID:
                continue
            MPID.append(doc.data["material_id"])
            struct = doc.structure.to_ase_atoms()
            mp_id = doc.data["material_id"]
            volume = doc.data["volume"]
            crystal = str(doc.data["symmetry"]["crystal_system"])
            # run_types = doc.run_types
            db.write(
                struct,
                sample=struct.get_chemical_formula(),
                material_id=mp_id,
                volume=volume,
                crystal_system=crystal,
                opt=False,
            )
            msg = f"Added {mp_id} {crystal} {volume:.2f} to database "
            print(msg)
        print(marquee("The following keys were added \n", 80, "*"))
        print(indent("name", 4))
        print(indent("material_id", 4))
        print(indent("volume", 4))
        print(indent("crystal_system", 4))
        return db


def update_dict(dict, **kwargs):
    """
    Update the dictionary with the provided keyword arguments.
    """
    for key, value in kwargs.items():
        if key in dict.keys():
            dict[key] = value
    return dict


def prepare_calculation(
    dbname: str,
    xc: list[str] | str = "PBE",
    kspacing: float = 0.2,
    encut: int = 450,
    calc_inplace: bool = False,
    **kwargs,
):
    """
    Prepare the calculation input files for the bulk structures.
    Parameters
    ----------
    dbname : str
        Name of the database.
    xc : str, optional
        Exchange-correlation functional, by default "PBE".
    kspacing : float, optional
        K-point spacing, by default 0.2.
    encut : int, optional
        Energy cutoff, by default 520.
    """

    calc = Vasp(
        setups="recommended",
        istart=0,
        icharg=2,
        ncore=4,
        kspacing=kspacing,
        ispin=1,
        encut=encut,
        gamma=True,
        ismear=0,
        sigma=0.05,
        lreal="Auto",
        isif=3,
        ibrion=2,
        nsw=500,
        algo="Veryfast",
        ediffg=-0.02,
        ediff=1e-5,
        prec="Normal",
        nelmin=4,
        nelm=120,
    )
    if isinstance(dbname, str):
        db = connect(dbname)
    else:
        db = dbname
    for row in db.select():
        params = row.key_value_pairs
        atoms = db.get_atoms(row.id)

        update_params = update_dict(params, **kwargs)
        if isinstance(xc, str):
            xc = [xc]
        if len(xc) == 0:
            raise ValueError("Please provide at least one xc functional.")

        for ppp in xc:
            atoms.calc = calc
            atoms.calc.set_xc_params(ppp)
            atoms.calc.set(gga=None)
            atoms.calc.set(
                directory="%s_%s_%s_%s"
                % (
                    update_params.get("sample", "unknown"),
                    update_params.get("crystal", "p1"),
                    update_params.get("xc", ppp),
                    update_params.get("material_id", "mp0000"),
                )
            )
            if calc_inplace:
                print(marquee(f" Running calculation for {row.material_id} ", 80, "*"))
                atoms.get_potential_energy()
                db.write(
                    atoms,
                    opt=True,
                    sp_energy=atoms.get_potential_energy(),
                    xc=ppp,
                    **update_params,
                )
            else:
                if not os.path.isdir(atoms.calc.directory):
                    try:
                        os.makedirs(atoms.calc.directory)
                    except FileExistsError as e:
                        msg = f"Directory {atoms.calc.directory} already exists."
                        raise RuntimeError(msg) from e

                atoms.calc.write_input(atoms)
                print(
                    marquee(
                        f" Calculation input files for {row.material_id} prepared ",
                        80,
                        "*",
                    )
                )


def check_converage(outcar_path):
    converage = False
    with open(outcar_path) as fd:
        lines = fd.readlines()
        for line in lines:
            if "reached required" in line:
                return True
        return converage


def get_energy_fromOszicar(oszicar_path):
    with open(oszicar_path, "r") as fd:
        lines = fd.readlines()
    return lines[-1].split()[4]


def update_db(dbname, workdir):
    if isinstance(dbname, str):
        db = connect(dbname)
    else:
        db = dbname
    paths = glob.glob(os.path.join(workdir, "**", "OUTCAR"), recursive=True)
    for path in paths:
        if not check_converage(path):
            print(f"{path.replace("OUTCAR", "")} not converged")
            continue
        fname = path.split("/")[-2].split("_")
        sample = fname[0]
        crystal = fname[1]
        xc = fname[2]
        material_id = fname[3]
        try:
            energy = np.loadtxt(os.path.join(path.replace("OUTCAR", "out.energy")))[0]
        except FileNotFoundError:
            energy = get_energy_fromOszicar(path.replace("OUTCAR", "OSZICAR"))
        atoms = read(path.replace("OUTCAR", "CONTCAR"))
        db.write(
            atoms,
            sample=sample,
            sp_energy=float(energy),
            material_id=material_id,
            crystal=crystal,
            xc=xc,
            opt=True,
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build the bulk database from materials project."
    )
    parser.add_argument(
        "--elements",
        type=list,
        required=True,
        help="List of chemical systems to query.",
    )
    parser.add_argument(
        "--database",
        type=str,
        required=True,
        help="Name of the database.",
    )
    parser.add_argument(
        "--xc",
        type=str,
        default="PBE",
        help="Exchange-correlation functional, by default 'PBE'.",
    )
    parser.add_argument(
        "--calc_inplace",
        type=bool,
        default=False,
        help="Whether to run the calculation in place.",
    )
    parser.add_argument(
        "--update_db",
        type=bool,
        default=False,
        help="Whether to update the database with new entries.",
    )

    parser.add_argument(
        "--workdir",
        type=str,
        default=".",
        help="Working directory for the calculations, by default current directory.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.elements is None or args.database is None:
        raise ValueError("Please provide the elements and database name.")

    docs = query_structure(
        chemsystem=[args.elements],
        fields=["volume", "material_id", "symmetry"],
    )
    db = parse_bulk_data(args.database, docs, args.xc)

    if not args.update_db:
        prepare_calculation(db, xc=args.xc, calc_inplace=args.calc_inplace)
    else:
        update_db(db, args.workdir)


if __name__ == "__main__":
    main()
