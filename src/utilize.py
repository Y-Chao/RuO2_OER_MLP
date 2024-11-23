#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'


import os
import sys
sys.path.insert(0, '/Users/ychao/Documents/develops/ase/')

import numpy as np
from ase import Atom, Atoms
from ase.io import read, write
from ase.db import connect
from ase.geometry import get_layers, get_distances
from ase.build import surface
from ase.build.surfaces_with_termination import surfaces_with_termination
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.data import covalent_radii, atomic_numbers
from ase.dev.nsurface import align_zdir


def rand_direction():
    """Generate a random vector, with distribution around the sphere
    Marsaglia method (https://mathworld.wolfram.com/SpherePointPicking.html)
    That is, picking (x1, x2) from uniform distribution on [-1, 1] respectivity and then
    rejects points for which x1^2 + x2^2 >= 1. From the remaining points,
        x = 2x_1 sqrt(1-x_1^2-x_2^2)
        y = 2x_2 sqrt(1-x_1^2-x_2^2)
        z = 1 - 2(x_1^2 + x_2^2)
    
    Returns:
    --------
    vector: array-like, shape (3,)
        A random vector
    """
    x1, x2 = 1, 1
    while x1**2 + x2**2 >= 1:
        x1 = np.random.uniform(-1, 1)
        x2 = np.random.uniform(-1, 1)
    return np.array(
        [
            2 * x1 * np.sqrt(1 - x1**2 - x2**2),
            2 * x2 * np.sqrt(1 - x1**2 - x2**2),
            1 - 2 * (x1**2 + x2**2),
        ]
    )

def calc_angle(vector, ref_vector = [0, 0, 1]):
    """
    Calculate the angle between the vector and the reference vector, the default
    reference vector is z-axis([0, 0, 1]).
    Parameters:
    -----------
    vector: array-like, shape (3,)
        The vector to be calculated
    ref_vector: array-like, shape (3,)
        The reference vector
    Returns:
    --------
    theta: float
        The angle (in degree) between the vector and the reference.
    """
    theta = np.arccos(np.dot(vector, ref_vector) / 
                     (np.linalg.norm(vector) * np.linalg.norm(ref_vector)))
    #print(vector, ref_vector, theta)
    return np.degrees(theta)

def constrain_direction(ref_vector, max_angle):
    """
    Generate a random vector within the a certain angle constraint.
    """
    v = rand_direction()
    iter = 0
    #print(ref_vector)
    while calc_angle(v, ref_vector) > max_angle:
        if iter > 100:
            raise ValueError('Fail to generate a random vector within the angle constraint!')
        v = rand_direction()
        iter += 1
    return v

def BLDA(elem_1, elem_2, sigma=0.1, scale=1):
    """
    Bond length distribution analysis (BLDA) is a method to generate the bond length
    distribution of a given element pair. The method is based on the Gaussian
    distribution of the bond length. The mean value of the Gaussian distribution is
    the covalent radius of the element pair, and the standard deviation is the
    parameter sigma. The default value of sigma is 0.1. The scale parameter is used
    to scale the bond length. ref: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00994
    """
    r = covalent_radii[elem_1] + covalent_radii[elem_2]
    return np.random.normal(r, sigma, 1)[0] * scale

def BLDA_vector(elem_1, elem_2, ref_vector=[0,0,1], max_angle=45):
    """
    Grow an atom with bond length distribution analysis (BLDA) method
    """
    vector = constrain_direction(ref_vector, max_angle)
    bond_length = BLDA(elem_1, elem_2)
    return bond_length * vector

def grow_atom_BLDA(surf, surf_atom_index=None, add_elem='H', ref_vector=[0,0,1], max_angle=45):
    """
    Grow an atom with bond length distribution analysis (BLDA) method
    """
    atoms = surf.copy()
    add_pos = []
    if surf_atom_index is None:
        raise ValueError('Please provide the index list of the surface atoms!')
    n_grow_atom = len(surf_atom_index)
    for i, index in enumerate(surf_atom_index):
        pos_0 = atoms.positions[index]
        accept = False
        iter = 0
        while iter < 100:
            pos_1 = pos_0 + BLDA_vector(atomic_numbers[atoms[index].symbol], atomic_numbers[add_elem], ref_vector, max_angle)
            _, dis = get_distances(pos_1, atoms.positions, cell=atoms.cell, pbc=True)
            iter += 1
            if np.all(dis > 0.5):
                accept = True
                break
        if accept:
            atoms.append(Atom(add_elem, pos_1))
            add_pos.append(pos_1)     
        else:
            print(f'Fail to grow atom {i} after {iter} iterations!') 
    return atoms, add_pos

def fix_layer(surf, fix_layer, miller=(0, 0, 1)):
    layer_index, dis_list = get_layers(surf, miller)
    n_layer = len(set(layer_index))
    c = FixAtoms(indices=np.argwhere(layer_index < fix_layer).flatten())
    tmp = surf.copy()
    tmp.set_constraint(c)
    return tmp

def parse_atom(atoms, *kwargs):
    """
    Parse one atom object to a dictionary, based on the given keywords.
    """
    results = {}
    if kwargs == None:
        return {}
    for k, v in atoms.__dict__.items():
        if k in kwargs:
            results.update({k: v})
    return results

def to_database(dbname, atoms, *kwargs):
    if len(atoms) == 0:
        atoms = [atoms]
    for atom in atoms:
        results = parse_atom(atom, *kwargs)
        with connect(dbname) as db:
            try:
                ene = atom.get_potential_energy()
            except:
                ene = None
            optimized = False
            db.write(atom, energy=ene, optim_flag=optimized , **results)
    return db

def extract_md(file, start_index=0, end_index=-1, num_structure=None, dbname=None):
    """
    Extract the atoms from the trajectory file.
    Parameters:
    -----------
    file: str
        The trajectory file.
    start_index: int
        The start index of the trajectory.
    end_index: int
        The end index of the trajectory.
    num_structure: int
        The number of structures to be extracted.
    dbname: str
        The database name.
    Returns:
    --------
    extract_traj: list
        The extracted atoms list.
    To do:
    ------
    1. Add the trajectory type judgement.
    2. Add the energy and force read.

    """
    trajectory = read(file, index='%s:%s'%(start_index, end_index))
    extract_traj = []
    if num_structure is None:
        step = int(len(trajectory)/200)
    for i, atoms in enumerate(trajectory):
        if i % step == 0:
            print(f'Extracting {i}th structure...')
            extract_traj.append(atoms)
    if dbname is not None:
        db = connect(dbname)
        for atoms in extract_traj:
            db.write(atoms)
    return extract_traj
