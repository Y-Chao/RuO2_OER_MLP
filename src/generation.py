#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'

"""
Define the generation section.
"""
import os
import re
from abc import ABC, abstractmethod
from uuid import uuid4
from datetime import datetime
import numpy as np
from ase.io import read, write
from ase.data import covalent_radii
from ase.geometry import get_distances
from candidates import BaseCandidate, Candidate

class Generation(ABC):

    default_task_type = {"growth": False, 
                 "dissolution": False, 
                 "reconstruction": True}
    default_limit = {"c1": 0.8, "c2": 1.2}
    
    def __init__(self, candidate=None, sources=None, dimension=3, c1=None, c2=None, use_mic=True, growth=False, dissolution=False, reconstruction=False, verbose=True, seed=9999):
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
        self.dimension = dimension
        self.set_limit(c1, c2)
        self.seed = seed # set the seed to generate the random number
        self.workdir = os.getcwd()
        self.verbose = verbose

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

    def _get_candidates(self, parent, indices=None):
        """
        Get the candidates.
        """
        return []
    
    def get_candidates(self, parent=None, indices=None):
        """
        Method used to get new candidates.
        """
        if parent is None:
            parent = self.template
        # if mol is None:
        #     mol = self.source_library

        if not parent.has_metadata('parent_search_iterations'):
            parent.add_metadata("parent_search_iterations", 0)
        if not parent.has_metadata('uuid'):
            parent.add_metadata("uuid", str(uuid4()))
 
        # generate new candidates
        candidates = self._get_candidates(parent, indices)
        for candiate in candidates:
            candiate.add_metadata("generator", self.name)
            candiate.add_metadata("uuid", str(uuid4()))
            candiate.add_metadata("parent_search_iterations", parent.get_metadata("parent_search_iterations") + 1)
            candiate.add_metadata("parent_uuid", parent.get_metadata("uuid"))
        
        if self.verbose:
            self.writer(candidates)
        
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

    def check_new_positions(self, candidate, new_positions, move_atomic_numbers, skip_index=None):
        """
        Check the new positions.
        Parameters:
        -----------
        candidate: Candidate
            The candidate object.
        new_positions: np.array
            The new positions of one atom.
        move_atomic_numbers: int
            The atomic number of the corresponding atom.
        skip_index: list
            The index list of the atom that should be skipped.
        """
        covalent_radii_move = covalent_radii[move_atomic_numbers]
        success = False
        for i in range(len(candidate)):
            if i in skip_index:
                continue
            covalent_radii_i = covalent_radii[candidate[i].number]
            covalent_bond_ij = covalent_radii_move + covalent_radii_i
            r_min = self.c1 * covalent_bond_ij
            r_max = self.c2 * covalent_bond_ij
            if self.mic:
                _, distances = get_distances(new_positions, candidate.positions[i], cell=candidate.cell, pbc=candidate.pbc)           
            else:
                distances = np.linalg.norm(new_positions - candidate.positions[i])
            if distances < r_min:
                return False
            elif not distances > r_max:
                success = True
        return success
    
    def writer(self, candidates):
        """
        Write the log file to a log_file, which named as the data_time.log
        """
        max_number = 0
        today = datetime.now()
        logname = today.strftime("%Y%m%d")  
        for file in os.listdir(self.workdir):
            match = re.match(r'\d{8}_(\d+).log', file)
            if match:
                number = int(match.group(1))
                max_number = max(max_number, number)

        logname = logname + "_%d.log"%(max_number + 1)
        with open(logname, 'w') as f:
            f.write("="*50)
            f.write('\n')
            f.write("Generation Log")
            f.write('\n')
            f.write("="*50)
            f.write('\n')
            f.write("The generation is %s \n"%self.name)
            f.write("The generations is started at %s\n"%today.strftime("%Y-%m-%d %H:%M:%S"))
            f.write("The seed is %d\n"%self.seed)
            f.write('\n')
            # write the header
            f.write(f"{'#Candidate UUID':^50}{'Parent UUID':^50}{'Parent Search Iterations':^30}{'Generator':^40}")
            f.write('\n')
            for candidate in candidates:
                f.write(f"{candidate.get_metadata('uuid'):^50}{candidate.get_metadata('parent_uuid'):^50}{candidate.get_metadata('parent_search_iterations'):^30}{candidate.get_metadata('generator'):^40}")
                f.write('\n')
            f.write("="*50)
            f.write('\n')
            f.write("The generations is ended at %s"%datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            f.write('\n')
    
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
        
        if  isinstance(mol, list) and len(mol) == 1:
            mol = [mol]
        
        for m in mol:
            for p in parent:
                pass
                
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
    
    def check_new_positions(self, candidate, new_positions, mol):
        """
        Check the new positions.
        """
        if isinstance(mol, list) and len(mol) == 1:
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
        

class Reconstruction_Generation(Generation):
    
    name = "ReconstructionGeneration"
    
    def __init__(self, atoms=None, rattle_amplitude=3, attempts=100, reconstruction=True, **kwargs):
        super().__init__(atoms, **kwargs)
        self.rattle_amplitude = rattle_amplitude
        self.attempts = attempts

    def get_indices_to_rattle(self, parent, indices=None):
        if indices is None:
            rattle_indices = parent.get_surface()
        else:
            rattle_indices = indices
        self.template.metadata['rattle_indices'] = rattle_indices
        return rattle_indices
    
    def _get_candidates(self, parent, indices=None):
        """
        Get the candidate by the reconstruction method.
        That is, the atom will be rattled by a random vector accordingly.
        """
        rattle_indices = self.get_indices_to_rattle(parent, indices)
        #np.random.seed(self.seed)
        candidate = parent.copy()
        accept_rattle_indices = []

        for i in rattle_indices:
            for _ in range(self.attempts):
                rattle_distance = np.random.normal(loc=0, scale=self.rattle_amplitude)
                rattle_vector = self.get_vector(rattle_distance)
                new_positions = parent.positions[i] + rattle_vector
                if self.check_new_positions(candidate, new_positions, parent.numbers[i], skip_index=[i]):
                    candidate[i].position = new_positions
                    accept_rattle_indices.append(i)
                    break
        candidate.add_metadata('actual_rattle_index', accept_rattle_indices)
        return [candidate]


class Dissolution_Generation(Generation):
    
    name = "DisslutioGeneration"
    
    def __init__(self, atoms=None, attempts=100, dissolution=True, **kwargs):
        super().__init__(atoms, **kwargs)
        self.attempts = attempts

    def _get_candidates(self, parent, indices=None):
        """
        Get the candidate by the dissolution method.
        """
        print(self.__dict__)
        dis_indices = self.get_indices_to_dis(parent, indices)
        candidate = parent.copy()
        for i in dis_indices:
            for _ in range(self.attempts):
                candidate.pop(i)
                new_bm = candidate.get_bondmatrix()
                conn = np.sum(new_bm, axis=0)
                if np.all(conn > 0):
                    break
        return [candidate]
    
    def get_indices_to_dis(self, parent, indices=None, num=1):
        if indices is None:
            dis_indices = parent.get_surface()
            # We will add some sorted method to get the one atom to be dissolved.
        else:
            dis_indices = indices
        return self.score_dissolution_index(parent, dis_indices, num)
        
    def score_dissolution_index(self, parent, dis_indices, num=1):
        """
        Score the dissolution based on the coordination number, z position, atomic number.
        First thing, we construct the coordination matrix.
        """
        dis_indices_coord = parent.coord_matrix[dis_indices]
        dis_numbers = parent.numbers[dis_indices]

        score_list = np.ones(len(dis_indices))
        for i in range(len(dis_indices)):
            # limit the metal atom and the coordinate saturation atom.
            if dis_numbers[i] < 18:
                score_list[i] += 1
            elif dis_indices_coord[i] == dis_indices_coord.min():
                score_list[i] += 1
        
        # sort the score list
        success_rate = np.array(score_list) / np.sum(score_list)
        dis_index = np.random.choice(dis_indices, size=num, p=success_rate)
        return dis_index
    
class Growth_Generation(Generation):

    name = "GrowthGeneration"
    

def main():
    path = '/Users/ychao/Library/CloudStorage/OneDrive-个人/nus/project/1_RuO2_stability/src/test/dateset/'
    candidate = Candidate.from_atoms(path +'RuO2_110_2x2_4L.vasp')
    candidate.get_surface()
    reconstruct = Reconstruction_Generation(candidate, attempts=100)
    new_candidate = reconstruct.get_candidates()
    if len(new_candidate) > 0:
        write(path + 'reconstruct.vasp', new_candidate)
    dissolution = Dissolution_Generation(candidate, attempts=100)
    new_candidate = dissolution.get_candidates()
    if len(new_candidate) > 0:
        write(path + 'dissolution.vasp', new_candidate)
    new_1 = Dissolution_Generation(new_candidate[0], attempts=100)
    new_candidate = new_1.get_candidates()
    if len(new_candidate) > 0:
        write(path + 'dissolution_1.vasp', new_candidate)


if __name__ == '__main__':
    main()

