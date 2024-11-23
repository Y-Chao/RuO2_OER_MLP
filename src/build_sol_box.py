#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'

import sys
import os
import argparse
import numpy as np
from ase.io import read, write
from ase.build import molecule
from ase.units import mol


def calc_sol_num(slab, sol, density=1.0, surf_height=2.0, scale=0.9):
    """Calculate the number of solvent molecules according to the density of the solvent.
    Parameters:
    -----------
    slab: ase.Atoms
        The slab surface
    sol: ase.Atoms
        The solution molecular
    density: float
        The density of the solvent, default is 1.0 g/mL
    surf_height: float
        The height of above the slab, default is 2.0 Ang
    Returns:
    --------
    num: int
        The number of solvent molecules
    """
    a, b, c = slab.cell
    assert np.allclose(slab.cell.cellpar()[3:5], 90.0), 'The slab is not orthogonal.'

    slab_area = np.linalg.norm(np.cross(a, b))  # units: Angstrom^2
    z = np.linalg.norm(c) - surf_height # units: Angstrom
    sol_volume = slab_area * z * 10**-24 # units: mL
    sol_mass = np.sum(sol.get_masses())  # units: g/mol
    num = int(sol_volume * scale * density / sol_mass * mol)

    print('The number of solvent molecules is: {}'.format(num))
    return num
    
def generate_water_box(a, b, c, num, output_file='water_box.pdb'):
    input_file = 'water_box.in'
    with open(input_file, 'w') as f:
        f.write('# The input file for solvation box build by packmol\n')
        f.write('tolerance 2.0\n')
        f.write('output {}\n'.format(output_file))
        f.write('structure sol.pdb\n')
        f.write('  number {}\n'.format(num))
        f.write('  inside box 0. 0. 0.  {} {} {}\n'.format(a, b, c))
        f.write('end structure\n')

    exit_status = os.system('packmol < {}'.format(input_file))
    if exit_status != 0:
        print('Error in conducting packmol. Please check the packmol installation.')
        for f in ['sol.pdb', 'water_box.pdb', 'water_box.in']:
            if os.path.exists(f):
                os.remove(f)
        return None
    
    water_box = read('water_box.pdb')
    for f in ['sol.pdb', 'water_box.pdb', 'water_box.in']:
        if os.path.exists(f):
            os.remove(f)
    return water_box

def add_water_box(slab, water_box, surf_height):
    """Add the water box above the slab surface.
    Parameters:
    -----------
    slab: ase.Atoms
        The slab surface
    water_box: ase.Atoms
        The water box
    Returns:
    --------
    slab_water: ase.Atoms
        The slab with water box
    """
    z_max = slab.positions[:, 2].max()
    water_box.positions[:, 2] = z_max + surf_height
    slab_water = slab + water_box
    return slab_water.wrap()

def main():
    parse = argparse.ArgumentParser(description='Build a solvation box above the slab surface.')
    parse.add_argument('slab_file', type=str, help='The file of slab surface.')
    parse.add_argument('-s', '--sol', type=str, help='The solvent to be used.')
    parse.add_argument('-n', '--num', type=int, help='The number of solvent molecules.')
    parse.add_argument('-w', '--water_box', type=str, help='The structure of water box.')
    parse.add_argument('-o', '--output', type=str, default='surf_water.vasp', help='The output file name.')

    args = parse.parse_args()

    slab = read(args.slab_file)
    a, b, c = slab.cell.cellpar()[:3]

    if args.water_box is not None:
        water_box = read(args.water_box)
    else:
        if args.sol is None:
            sol = molecule('H2O')
            write('sol.pdb', sol)
        else:
            sol = read(args.sol)
    
        if args.num is None:
            num = calc_sol_num(slab, sol, density=1.0)
        else:
            num = args.num

        water_box = generate_water_box(a, b, c, num)
    
    if water_box is None:
        print('Fail to generate the water box!')
        sys.exit(1)
    
    atoms = add_water_box(slab, water_box, surf_height=2.0)
    write(args.output, atoms)


if __name__ == '__main__':
    main()
