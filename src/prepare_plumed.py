#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = 'Chao Yang'
__version__=	'1.0'

import sys
from ase.io import read, write

def main():
    # read the strucutre file
    atoms = read(sys.argv[1])

    # define the region
    region_1 = sys.argv[2]
    region_2 = sys.argv[3]
    print("The surface atoms is belowe the region: ", region_1)
    print("The water atoms is above the region: ", region_2)

    print()


if __name__ == '__main__':
    main()