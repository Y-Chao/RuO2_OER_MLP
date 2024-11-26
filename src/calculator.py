#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'


"""
Here is a simple calculator module used to calculate the energy, forces or magmom of a structure.
Meanwhile, it can also be interfaced with the optimization and md modules.
Two main calculators are defined here: VASP and CHGNet.
"""

import warnings
from ase import Atoms
from ase.calculators.vasp import Vasp
from chgnet.model.model import CHGNet
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.io.ase import AseAtomsAdaptor
from chgnet.model.dynamics import MolecularDynamics
from chgnet.model import StructOptimizer


class BaseCalculator:
    def __init__(self, calc_name='VASP', task_type='sp', **kwargs) -> None:
        self.calc_name = calc_name
        self.task_type = task_type
        self.parameters = kwargs

    def initialize_calculator(self):
        if self.calc_name == 'VASP':
            return Vasp()
        elif self.calc_name == 'CHGNet':
            return CHGNet.load()
    
    def calculate(self, atoms):
        pass

class VaspCalculator(BaseCalculator):
    def __init__(self, **kwargs) -> None:
        super().__init__(calc_name='VASP', **kwargs)
        self.calc = self.initialize_calculator()
    
    def set_calculatr_parameters(self, atoms, **kwargs):
        self.calc.set(**kwargs)
        atoms.set_calculator(self.calc)
        return atoms

    def run(self, atoms):
        if self.task_type == 'sp':
            energy = atoms.get_potential_energy()
            return atoms.positions, energy
        elif self.task_type == 'opt':
            optimizer, fmax, opt_steps = self.parameters['optimizer'], self.parameters['fmax'], self.parameters['opt_steps']
            dyn = optimizer(atoms)
            dyn.run(fmax=fmax, steps=opt_steps)
            return atoms.positions, atoms.get_potential_energy()
        
    def set_optimizer(self):
        optimizer = self.parameters.get('optimizer', 'BFGS')
        fmax = self.parameters.get('fmax', 0.02)
        steps = self.parameters.get('steps', 500)
        return optimizer, fmax, steps

class ChgnetCalculator(BaseCalculator):
    def __init__(self, **kwargs) -> None:
        super().__init__(calc_name='CHGNet', **kwargs)
        self.calc = self.initialize_calculator()

    def atoms2structure(self, atoms):
        if isinstance(structure, Structure):
            structure = structure
        elif isinstance(structure, Atoms):
            structure = AseAtomsAdaptor.get_structure(structure)
    
    def run(self, atoms):
        atoms = self.atoms2structure(atoms)
        if self.task_type == 'sp':
            prediction = self.calc.predict_structure(atoms)
            for key, unit in [
                ("energy", "eV/atom"),
                ("forces", "eV/A"),
                ("stress", "GPa"),
                ("magmom", "mu_B"),
                ]:
                print(f"CHGNet-predicted {key} ({unit}):\n{prediction[key[0]]}\n")
            return atoms.positions, prediction['energy']
        
        elif self.task_type == 'opt':
            relaxer = StructOptimizer()
            result = relaxer.relax(atoms)
            print("CHGNet relaxed structure", result["final_structure"])
            print("relaxed total energy in eV:", result['trajectory'].energies[-1])
            return result['trajectory'].energies[-1]
        
        elif self.task_type == 'md':
            md = MolecularDynamics(
                atoms=atoms,
                model=self.calc,
                ensemble=self.parameters['ensemble'],
                temperature=self.parameters['temperature'],  # in K
                timestep=self.parameters['timestep'],  # in femto-seconds
                trajectory="md_out.traj",
                logfile="md_out.log",
                loginterval=self.parameters['interval'],
                )
            md.run(self.parameters['steps'])
    

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




