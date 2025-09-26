#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

"""
Load toml parameters for building surface and interface models.
"""

import os
from argparse import Namespace

import tomllib


def load_surface_toml(file: str):
    """
    Load a toml files and check the configurations
    Parameters:
    -----------
    file: str
    """
    assert os.path.isfile(file), "The surface configuration file can not be found!"

    Default_bulk_parameters = ["bulk_db", "key_info"]
    Default_surface_parameters = [
        "terminations",
        "miller_index",
        "layers",
        "vacuum",
        "min_lattice",
    ]
    Default_solvation_parameters = [
        "solvation",
        "pH",
        "ions",
        "ions_number",
    ]

    with open(file, "br") as fd:
        config = tomllib.load(fd)

    for key, value in config.items():
        if key == "bulk_info":
            missing_key = set(Default_bulk_parameters) - set(value.keys())
            if missing_key:
                raise KeyError(
                    f"bulk_info is missing the following keys: {missing_key}"
                )
        elif key == "surface_info":
            missing_key = set(Default_surface_parameters) - set(value.keys())
            if missing_key:
                raise KeyError(f"surface is missing the following keys: {missing_key}")

            if len(value["terminations"]) == 0 and isinstance(
                value["terminations"], list
            ):
                value["terminations"] = None
            if len(value["miller_index"]) != 3:
                raise ValueError("the miller index needs three integral.")
            if value["min_lattice"] < 6:
                raise ValueError(
                    f"to keep the surface area enough, the min_lattice {value["min_lattice"]} is too small."
                )
        elif key == "ion_info":
            missing_key = set(Default_solvation_parameters) - set(value.keys())
            if missing_key:
                raise KeyError(f"surface is missing the following keys: {missing_key}")
            if len(value["ions"]) == 0 and isinstance(value["ions"], list):
                value["ions"] = None
            if len(value["ions_number"]) == 0 and isinstance(
                value["ions_number"], list
            ):
                value["ions_number"] = None

            if len(value["ions"]) != len(value["ions_number"]):
                raise KeyError("the ions number should equals to ions")
        else:
            pass
    return Namespace(**config)


def main():
    args = load_surface_toml(
        "/Users/c_yang/Library/CloudStorage/OneDrive-Personal/nus/project/1_RuO2_stability/examples/surface_configuration.toml"
    )


if __name__ == "__main__":
    main()
