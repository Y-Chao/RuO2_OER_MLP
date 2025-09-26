#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

from dataclasses import dataclass


@dataclass
class Cp2kInfo:
    version: str = None
    restart: bool = None
    terminated_by_request: bool = None


CP2K_INFO_VERSION_PATTERN = r"""(?xm)
    

    """
