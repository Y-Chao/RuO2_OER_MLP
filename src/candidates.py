#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'

"""
Define the candidate class, which enherits from the Atoms class.

"""

import os
from abc import ABC
from copy import deepcopy
from ase import Atoms
from ase.io import read, write

class BaseCandidate(ABC, Atoms):

    def __init__(self, atoms):
        """
        Class for candidate, which enherits from the Atoms class.
        Parameters:
        atoms: Atoms
            The atoms object.
    
        """
        kwargs = self.parse_atoms(atoms)
        Atoms.__init__(self, **kwargs)

    def parse_atoms(self, atoms):
        """
        Parse the atoms object to get the information.
        """
        kwargs = atoms.todict()
        if 'constraints' in kwargs:
            constraints = kwargs.pop('constraints')
            kwargs['constraint'] = constraints
        return kwargs
    
    def compare(self, identifier):
        for a, b in zip(identifier, self.get_identifier()):
            equal = ( a == b ).all
            if not equal:
                return equal
        return equal
    
    def get_identifier(self):
        return (self.get_atomic_numbers(), self.get_positions(), self.get_cell(), self.get_pbc())

    def copy(self):
        """
        Copy the candidate.
        """
        kwarg = self.to_dict().copy()
        metadata = self.metadata.copy()
        new_candidate = self.__class__(**kwarg)
        new_candidate.metadata = metadata
        new_candidate.set_constriants(deepcopy(self.get_constraints()))
        return new_candidate
    
    def has_metadata(self, key):
        """
        Check if the metadata has the key.
        Parameters:
        -----------
        key: str
            The key of the metadata.
        """
        return key in self.metadata
    
    def add_metadata(self, key, value):
        """
        Add an entry to metadata
        Parameters:
        -----------
        key: str
            The key of the metadata.
        value: any
            The value of the metadata.
        """
        self.metadata[key] = value

    def get_metadata(self, key):
        """
        Get the metadata by key.
        Parameters:
        -----------
        key: str
            The key of the metadata.
        """
        return self.metadata.get(key, None)
    
    def pop_metadata(self, key):
        """
        Pop the metadata by key.
        Parameters:
        -----------
        key: str
            The key of the metadata.
        """
        return self.metadata.pop(key, None)
    
    def reset_metadata(self):
        """
        Reset the metadata.
        """
        self.metadata = dict()

    def set_calculator(self, calculator):
        """
        Set the calculator.
        Parameters:
        -----------
        calculator: Calculator
            The calculator object.
        """
        self.calc = calculator

class Candidate(BaseCandidate):
    """
    Class for candidate, which enherits from the Atoms class.
    """

    @classmethod
    def from_atoms(self, atoms):
        """
        Class for candidate, which enherits from the Atoms class.
        Parameters:
        atoms: Atoms
            The atoms object.
    
        """
        if isinstance(atoms, str):
            if os.path.exists(atoms):
                atoms = read(atoms)
            else:
                raise ValueError("The file %s does not exist."%atoms)
        else:
            atoms = atoms
        return self(atoms)
            

def main():
    path = '/Users/ychao/Library/CloudStorage/OneDrive-个人/nus/project/1_RuO2_stability/src/test/dateset/RuO2_110_2x2_4L.vasp'
    candidate = Candidate.from_atoms(path)
    print(candidate)

if __name__ == '__main__':
    main()
