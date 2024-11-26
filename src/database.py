#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'


from abc import ABC
import sqlite3

def BaseDB(ABC):
    def __init__(self, dbname):
        self.dbname = dbname
    
    def write(self, *args, **kwargs):
        pass

class DataBase(BaseDB):
    initial_state = [
    ]

    def __init__(self, dbname=None, initialize=True):
        super().__init__(dbname)
        if initialize and os.path.exists(self.dbname):
            os.remove(self.dbname)
        
        self.connect = sqlite3.connect(self.dbname,  timeout=600)
        self.write_interval = write_interval

        self.storage_keys = [
            'positions',
            'energy',
            'forces',
            'cell',
            'pbc',
            'iteration',
            'types',
            'uuid'
        ]

        
