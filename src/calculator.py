#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'


"""
Here is a simple calculator module used to calculate the energy, forces or magmom of a structure.
Meanwhile, it can also be interfaced with the optimization and md modules.
"""

import warnings
from ase import Atoms
from chgnet.model.model import CHGNet
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.io.ase import AseAtomsAdaptor
from chgnet.model.dynamics import MolecularDynamics
from chgnet.model import StructOptimizer


def static_calculator(structure):
    """
    Calculate the energy, forces and stress of a structure by CHGNet.
    """
    if isinstance(structure, Structure):
        structure = structure
    elif isinstance(structure, Atoms):
        structure = AseAtomsAdaptor.get_structure(structure)

    chgnet = CHGNet.load()
    prediction = chgnet.predict_structure(structure)

    for key, unit in [
        ("energy", "eV/atom"),
        ("forces", "eV/A"),
        ("stress", "GPa"),
        ("magmom", "mu_B"),
    ]:
        print(f"CHGNet-predicted {key} ({unit}):\n{prediction[key[0]]}\n")
    return prediction['energy']


def md_calculator(structure, temperature, timestep=1, steps=10000, interval=100, ensemble='nvt'):
    warnings.filterwarnings("ignore", module="pymatgen")
    warnings.filterwarnings("ignore", module="ase")

    if isinstance(structure, Structure):
        structure = structure
    elif isinstance(structure, Atoms):
        structure = AseAtomsAdaptor.get_structure(structure)

    
    chgnet = CHGNet.load()

    md = MolecularDynamics(
        atoms=structure,
        model=chgnet,
        ensemble=ensemble,
        temperature=temperature,  # in K
        timestep=timestep,  # in femto-seconds
        trajectory="md_out.traj",
        logfile="md_out.log",
        loginterval=interval,
        )
    md.run(steps)  # run a 0.1 ps MD simulation

def geo_calculator(structure):
    """
    Optimize the structure by CHGNet.
    """
    if isinstance(structure, Structure):
        structure = structure
    elif isinstance(structure, Atoms):
        structure = AseAtomsAdaptor.get_structure(structure)

    relaxer = StructOptimizer()
    result = relaxer.relax(structure)
    print("CHGNet relaxed structure", result["final_structure"])
    print("relaxed total energy in eV:", result['trajectory'].energies[-1])
    return result['trajectory'].energies[-1]




