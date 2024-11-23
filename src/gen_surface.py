#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'

"""
Although pymatgen surface provides great convenience to generate slabs, 
it may be not suitable for us to generate it.
"""

from pymatgen.core.surface import generate_all_slabs

def gen_slab_simple(bulk, miller_index, slab_layer, vacuum_layer):
    surface = generate_all_slabs(structure=bulk,
                                 miller_index=miller_index,
                                 min_slab_size=slab_layer,
                                 min_vacuum_size=vacuum_layer,
                                 max_normal_search=2,
                                 in_unit_planes=True
                                )
    return surface
    
