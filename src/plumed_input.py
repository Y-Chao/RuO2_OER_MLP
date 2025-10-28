#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"


"""
This module provides functions to generate PLUMED input files.

Generally, the plumed input file consists of four parts:
1. Global settings
    * Units:
        - Length: nm (default),  A (for Angstrom), um (for micrometer), Bohr(0.052917721067 nm)
        - Energy: kj/mol (default), j/mol, kcal/mol (4.184 kj/mol), eV (96.485 kj/mol), Hartree(2625.49962 kJ/mol)
        - Time: ps (default), fs, ns, atomic (2.418884326505e-05 ps)
        - Mass: amu (default)
        - Charge: e (default)
    * Group definitions:
        - o: GROUP ATOMS=1-10
            ...
        - p: GROUP ATOMS=11-20
            ...
2. Collective variable (CV) definitions
    s_co: COORINATION GROUPA=o GROUPB=p R_0=2.0 NN=12 NM=24 D_0=0.0
3. Biasing methods
    opes: OPES_METAD ARG=s_co PACE=500 SIGMA=0.1 HEIGHT=1.0 FILE=HILLS TEMP=300.0
4. Output settings
    PRINT ARG=s_co,opes.bias FILE=colvar STRIDE=100
"""

# group_dict = {} # A dictionary to store group definitions
# cv_dict = {}   # A dictionary to store collective variable definitions
# bias_dict = {} # A dictionary to store biasing method definitions
# print_dict = {} # A dictionary to store print definitions
