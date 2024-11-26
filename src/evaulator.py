#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'


from abc import ABC
from ase.optimize.bfgs import BFGS
from ase.optimize.fire import FIRE


class BaseEvaluator(ABC):
    def __init__(self, calculator, optimizer=BFGS, optimizer_kwargs={"fmax": 0.02, "steps":500}, store_trajectory=True):
        self.calculator = calculator
        self.optimizer = optimizer
        self.optimizer_kwargs = optimizer_kwargs
        self.store_trajectory = store_trajectory
