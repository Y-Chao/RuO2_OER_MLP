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
import ctypes
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry import get_layers

EXTERNAL_LIB = os.path.dirname(os.path.abspath(__file__))

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
        self.reset_metadata()

    def parse_atoms(self, atoms=None):
        """
        Parse the atoms object to get the information.
        """
        if atoms is None:
            atoms = self
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
        kwarg = self.todict().copy()
        metadata = self.metadata.copy()
        new_candidate = Candidate(self)
        try:
            new_candidate.metadata = metadata
        except AttributeError:
            new_candidate.reset_metadata()
        new_candidate.constriants = deepcopy(self.constraints)
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
    def from_atoms(cls, atoms):
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
        candidate = cls(atoms)
        return candidate
    
    def get_surface(self, atoms=None, miller=(0,0,1), deeper=False):
        """
        Get the surface of the candidate.
        """
        if atoms is None:
            atoms = self
        layer_index, _ = get_layers(atoms, miller)
        n_layer = len(set(layer_index))
        # To do: maybe the 1st layer is not enough to describe the surface.
        surface_index = [i for i, index in enumerate(layer_index) if index == n_layer-1]
        surface_metal = False
        try_layer = 0
        while not surface_metal:
            if np.all(atoms.numbers[surface_index] < 18):
                try_layer += 1
                self.add_metadata('surface_ads_index', deepcopy(surface_index))
                surface_index.extend([i for i, index in enumerate(layer_index) if index == n_layer-try_layer-1])
            else:
                surface_metal = True
        if deeper:
            surface_index.append([i for i, index in enumerate(layer_index) if index == n_layer-try_layer-2]) 
        try:
            self.add_metadata('surface_sub_index', list(set(surface_index) - set(self.metadata['surface_ads_index'])))
        except KeyError:
            self.metadata['surface_ads_index'] = []
            self.metadata['surface_sub_index'] = deepcopy(surface_index)
        self.add_metadata('surface_index', surface_index)
        return surface_index

    @property
    def bondmatrix(self):
        """
        Get the bond matrix of the candidate.
        """
        return self.get_bondmatrix()

    @property
    def distance_matrix(self):
        return self._distance_matrix

    def _distance_matrix(self):
        """
        Get the distance matrix of the candidate.
        """
        if np.sum(self.pbc) > 0:
            return self.get_distance_matrix(mic=True)
        else:
            return self.get_distance_matrix(mic=False)

    def calctypes(self):
        """
        Get the ctypes parameters for external lib.        
        """
        # conver the coordinates to the one-dimensional.
        coords_1d = self.positions.flatten()
        cell_1d = self.cell.flatten()
        frac_coords_1d = self.get_scaled_positions().flatten()
        bmatrix_size = len(self) ** 2
        
        c_natoms = ctypes.pointer(ctypes.c_int(len(self)))
        c_coords = ctypes.pointer((ctypes.c_double * len(coords_1d))(*coords_1d))
        c_cell = ctypes.pointer((ctypes.c_double *len(cell_1d))(*cell_1d))
        c_numbers = ctypes.pointer((ctypes.c_int *len(self))(*self.numbers))
        c_frac_coords = ctypes.pointer((ctypes.c_double * len(frac_coords_1d))(*frac_coords_1d))
        c_bondmatrix = ctypes.pointer((ctypes.c_int * bmatrix_size)(*[0 for i in range(bmatrix_size)]))
        return c_natoms, c_coords, c_cell, c_numbers, c_frac_coords, c_bondmatrix
    
    @property
    def coord_matrix(self):
        return self.get_coordination()

    def get_coordination(self):
        if not hasattr(self, 'bondmatrix'):
            bm = self.get_bondmatrix()
        else:
            bm = self.bondmatrix

        coord_matrix = np.zeros(len(self))
        for i in range(len(self)):
            coord_matrix[i] = np.sum(bm[i])
        
        return coord_matrix

    @property
    def bondneed(self):
        return self.get_bondneed()
    
    def get_bondneed(self):
        """
        Get the bond need of the candidate. This is not accurate for metal coordinates.
        We should further improve this function.
        """
        if not hasattr(self, 'coord_matrix'):
            coord_matrix = self.get_coordination()
        else:
            coord_matrix = self.coord_matrix
        satuated_bond = []
        for i in range(len(self)):
            satuated_bond.append(np.min([abs(self.numbers[i]-2), abs(self.numbers[i]-10), abs(self.numbers[i]-18)]))
        return np.array(satuated_bond) - coord_matrix

    def get_bondmatrix(self):
        """
        Get the bond matrix of the candidate.
        Here, we use the extern lib (fortran) to calculate the bond matrix.
        Return:
        -------
        bondmatrix: list
            The bond matrix.
        """
        external_lib = EXTERNAL_LIB + '/libs/checkminbond.so'
        
        c_natoms, c_coords, c_cell, c_numbers, c_frac_coords, c_bondmatrix = self.calctypes()
        bcal = ctypes.cdll.LoadLibrary(external_lib)
        bcal.get_all_bondmatrix_(c_natoms, c_numbers, c_coords, c_bondmatrix, c_cell)
        bm = np.array(c_bondmatrix.contents).reshape(len(self), len(self))
        return bm

def main():
    path = '/Users/ychao/Library/CloudStorage/OneDrive-个人/nus/project/1_RuO2_stability/src/test/dateset/RuO2_110_2x2_4L.vasp'
    candidate = Candidate.from_atoms(path)
    candidate.get_surface()
    print(candidate.coord_matrix[0], np.argwhere(candidate.bondmatrix[0]==1))
        

if __name__ == '__main__':
    main()
