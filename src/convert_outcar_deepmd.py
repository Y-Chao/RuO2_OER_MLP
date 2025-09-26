#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = 'Chao Yang'
__version__=	'1.0'

"""
Convert OUTCAR file to DeepMD format.

Usage:
    python convert_outcar_deepmd.py root_dir
"""

import os
import glob
import argparse
from tqdm import tqdm
import dpdata

def parse_args():
    parser = argparse.ArgumentParser(description="Convert OUTCAR file to DeepMD format.")
    parser.add_argument("-d", "--dir", type=str, help="Root directory containing OUTCAR files.")
    parser.add_argument("-s", "--suffix", type=str, default="OUTCAR", help="Suffix of the OUTCAR files to be processed.")
    return parser.parse_args()

def main():
    args = parse_args()
    root_dir = args.dir
    suffix = args.suffix
    paths = glob.glob(os.path.join(root_dir, f"*.{suffix}"), recursive=True)
    data = dpdata.MultiSystems()

    for path in tqdm(paths, total=len(paths), desc="Processing OUTCAR files"):    
        dsys = dpdata.LabeledSystem(path, fmt="vasp/outcar", type_map=['H', 'O', 'Ru'], begin=1)
        if dsys:
            data.append(dsys)
        data.to_deepmd_raw(f"{path.replace(f'.{suffix}', '')}")

if __name__ == '__main__':
    main()