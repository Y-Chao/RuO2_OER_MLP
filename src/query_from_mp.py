#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'

import json
from ase.db import connect
from mp_api.client import MPRester

HOME = '/Users/ychao'
def load_api_key(fname='%s/.mpg_api'%HOME):
    with open(fname, 'r') as f:
        jdata = json.load(f)
    return jdata['mp_api']


def query_structure(chemsystem, volume, fields=None):
    api_mpg = load_api_key()
    with MPRester(api_mpg) as mpr:
        docs = mpr.materials.search(
            chemsys=chemsystem,
            volume=volume,
            fields=fields
        )
    return docs

def parse_bulk_data(dbname, docs):
    with connect(dbname) as db:
        for doc in docs:
            struct = doc.structure.to_ase_atoms()
            mp_id = doc.material_id
            volume = doc.volume
            #run_types = doc.run_types
            db.write(struct, material_id=mp_id, volume=volume)
    return db
 
def main():
    docs = query_structure(
        chemsystem=['Ru-O', 'Ru-O-H', 'Ru'],
        volume=(0, 1000),
        fields=['structure', 'material_id', 'volume', 'task_ids', 'run_types']
        )
    parse_bulk_data('Ru_O_bulk.db', docs)

if __name__ == '__main__':
    main()