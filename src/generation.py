#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'

"""
Define the generation section.
"""
from abc import ABC, abstractmethod
from uuid import uuid4
import numpy as np
from ase.data import covalent_radii
from candidates import BaseCandidate

class Generation(ABC):

    default_task_type = {"growth": False, 
                 "dissolution": False, 
                 "reconstruction": True}
    default_limit = {"c1": 0.8, "c2": 1.2}
    
    def __init__(self, candidate=None, sources=None, use_mic=True, growth=False, dissolution=False, reconstruction=False, seed=9999):
        """
        Initialize the generation class.
        Parameters:
        -----------
        atoms: Atoms
            The atoms object. It act as the template to add, delete or manipulate the structure.
        source: list
            The source list. It defines the species that will be manipulated.
        use_mic: bool
            Whether to use minimum image convention.
        """
        self.template = candidate.copy()
        self.mic = use_mic
        self.source_library = []
        self.seed = seed # set the seed to generate the random number
        # define three types of generation
        task_type_number = int(growth + dissolution + reconstruction)
        if task_type_number != 1:
            self.task_type = self.default_task_type
        else:
            self.task_type = {"growth": growth,
                              "dissolution": dissolution,
                              "reconstruction": reconstruction
                              }
            
        if sources is None or len(sources) == 0:
            assert self.task_type['reconstruction'], "Only reconstruction can be conducted with no sources."
        else:
            assert not self.task_type['reconstruction'], "Only growth or dissolution can be conducted with sources."
            for mol in sources:
                is_exist = self.check_source(i)
                assert is_exist and self.task_type['dissolution'], "The source is not in initial atoms, it can not be dissolved."
                self.source_list.append(self.source_2_atoms(mol))
            
    @classmethod
    def source_2_atoms(mol):
        """
        Convert the source to atoms object.
        """
        from ase.build import molecule
        if isinstance(mol, Atoms):
            return mol
        elif isinstance(mol, str):
            return molecule(mol)
        else:
            raise ValueError("The %s type element in source is not supported."%type(mol))     

    def _get_candidates(self, parent):
        """
        Get the candidates.
        """
        return []
    
    def get_candidates(self, parent=None, mol=None):
        """
        Method used to get new candidates.
        """
        if parent is None:
            parent = self.template
        if mol is None:
            mol = self.source_library
        
        if not parent.has_metadata('parent_search_iterations'):
            parent.add_metadata("search_iterations", 0)
        if not parent.has_metadata('uuid'):
            parent.add_metadata("uuid", str(uuid4()))
 
        # generate new candidates
        candidates = self._get_candidates(parent, mol)
        for candiate in candidates:
            candiate.add_metadata("generator", self.name)
            candiate.add_metadata("uuid", str(uuid4()))
            candiate.add_metadata("parent_search_iterations", [p.get_metadata("search_iterations") for p in parent])
            candiate.add_metadata("parent_uuid", [p.get_metadata("uuid") for p in parent])
        return candidates
    
    @property
    @abstractmethod
    def name(self):
        pass

    def set_limit(self, c1=None, c2=None):
        """
        Set the limit of the bond length.
        """
        if c1 is None:
            c1 = self.default_limit['c1']
        if c2 is None:
            c2 = self.default_limit['c2']
        self.c1 = c1
        self.c2 = c2

    def check_new_positions(self, candidate, new_positions, mol):
        """
        Check the new positions.
        """
        if len(mol) == 1:
            mol = [mol]
        
        success = False
        for i in range(len(mol)):
            for j in len(candidate):
                covalent_radii_ij = covalent_radii[mol[i].number] + covalent_radii[candidate[j].number]
                r_min = self.c1 * covalent_radii_ij
                r_max = self.c2 * covalent_radii_ij
                if self.mic:
                    distances = candidate.get_mic(new_positions, candidate.positions[j], cell=candidate.cell, pbc=candidate.pbc)
                else:
                    distances = np.linalg.norm(new_positions - candidate.positions[j])
                
                if distances < r_min:
                    return False
                elif not distances > r_max: # which guranatees the new positions at least have one bond
                    success = True
        return success
        

class Growth_Generation(Generation):

    name = "GrowthGeneration"

    def __init__(self, atoms=None, source=None, dimension=3, use_mic=True):
        super().__init__(atoms, source, use_mic)
        self.dimension = dimension

    def _get_candidate(self, parent, mol):
        """
        Get the candidate by the growth method.
        Parameters:
        -----------
        Parents: Candidate
            The parent candidate.
        mol: Atoms
            The growth molecules.
        Returns:
        --------
        candidates: list
            The list of candidates.
        """
        candidates = []

    def get_sphere_vector(self, atomic_number_i, atomic_number_j):
        """
        Get the vecort length by the sum of two colvalent radius.
        Parameters:
        -----------
        atomic_number_i: int
            The atomic number of the first atom.
        atomic_number_j: int
            The atomic number of the second atom.
        Returns:
        --------
        vector: np.array
            The vector length.
        """ 
        dim = self.dimension
        covalent_bondlenght = covalent_radii[atomic_number_i] + covalent_radii[atomic_number_j]
        r_min = self.c1 * covalent_bondlenght
        r_max = self.c2 * covalent_bondlenght
        np.random.seed(self.seed)
        # the power of the r could make the uniform distribution
        r = np.random.uniform(r_min ** dim, r_max ** dim)**(1 / dim )
        return self.get_vector(r)
    
    def get_vector(self, r):
        """
        Baesed on the dimension limits, we get the sphere coordinate of the vector.
        And then, it was changed to the cartesian coordinate.
        Parameters:
        -----------
        r: float
            The length of the vector.
        Returns:
        --------
        vector: np.array
            The vector.

        To do:
        ------
        Update the sp_coord limitations for 1 and 2.
        """
        dim = self.dimension
        # spehere coordinate limits,
        #    1. Only on the X-positive axis
        #    2. Only on the XY plane
        #    3. On the whole sphere
        sp_coord = {
            1: {"theta": [0, 0], "phi":[np.pi/2, np.pi/2]},
            2: {"theta": [0, 2 * np.pi], "phi":[np.pi/2, np.pi/2]},
            3: {"theta": [0, 2 * np.pi], "phi":[0, np.pi]}
        }

        np.random.seed(self.seed)
        theta = np.random.uniform(*sp_coord[dim]["theta"])
        phi = np.random.uniform(*sp_coord[dim]["phi"])
        vector = r * np.array([np.sin(phi) * np.cos(theta),  # x = r sin(phi) cos(theta)
                               np.sin(phi) * np.sin(theta),  # y = r sin(phi) sin(theta)
                               np.cos(phi)]                  # z = r cos(phi)
                               )
        return vector
    
    

        


class Reconstruction_Generation(Generation):
    
    name = "ReconstructionGeneration"
    
    def __init__(self, atoms=None, source=None, use_mic=True):
        super().__init__(atoms, source, use_mic)
    
    def _get_candidate(self):
        pass


class Dissolution_Generation(Generation):
    
    name = "DisslutioGeneration"
    
    def __init__(self, atoms=None, source=None, use_mic=True):
        super().__init__(atoms, source, use_mic, growth=True)
    
    def _get_candidate(self):
        pass


