#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

"""
Load toml parameters for building surface and interface models.
"""

import os
import tomllib
from types import SimpleNamespace


def load_surface_toml(file: str):
    """
    Load a toml files and check the configurations
    Parameters:
    -----------
    file: str
    """
    assert os.path.isfile(file), "The surface configuration file can not be found!"

    Default_bulk_parameters = ["bulk_db", "crystal", "xc"]
    Default_surface_parameters = [
        "terminations",
        "miller_index",
        "layers",
        "vacuum",
        "min_lattice",
        "symmetry",
        "fix",
    ]
    Default_solvation_parameters = [
        "solvation",
        "num_sol",
        "surface_height",
        "pH",
        "ions",
        "ions_number",
        "verbose",
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
                if "layer_criteria" in missing_key:
                    config.get("layer_criteria", 0.5)
                elif "fix" in missing_key:
                    if value["symmetry"]:
                        config.get("fix", None)
                    else:
                        config.get("fix", "bottom")
                else:
                    raise KeyError(
                        f"surface is missing the following keys: {missing_key}"
                    )

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
        elif key == "interface_info":
            missing_key = set(Default_solvation_parameters) - set(value.keys())
            if missing_key:
                if "solvation" in missing_key:
                    config.get("solvation", None)
                elif "num_sol" in missing_key:
                    config.get("num_sol", None)
                elif "surface_height" in missing_key:
                    config.get("surface_height", 1.0)
                elif "pH" in missing_key:
                    config.get("pH", 7)
                elif "verbose" in missing_key:
                    config.get("verbose", False)
                elif "ions" in missing_key and "ions_number" in missing_key:
                    config.get("ions", None)
                    config.get("ions_number", None)
                else:
                    raise KeyError(
                        f"interface is missing the following keys: {missing_key}"
                    )
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
    return dict2ns(config)


def dict2ns(d):
    if isinstance(d, dict):
        return SimpleNamespace(**{k: dict2ns(v) for k, v in d.items()})
    elif isinstance(d, list):
        return [dict2ns(x) for x in d]
    else:
        return d


def main():
    args = load_surface_toml(
        "/Users/c_yang/Library/CloudStorage/OneDrive-Personal/nus/project/1_RuO2_stability/examples/xMOI_configurations.toml"
    )
    print(args)


if __name__ == "__main__":
    main()
