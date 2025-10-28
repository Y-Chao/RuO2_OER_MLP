#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

import sys
from pathlib import Path


class Cp2kOutput:
    """
    A class to parse the CP2K output file referred to the following link:
    https://robinzyb.github.io/cp2kdata/.

    Parameters:
    -----------
    file: str
        The CP2K output file.

    Returns:
    --------
    None
    """

    def __init__(self, file: str | Path, run_type: str = "md"):
        if isinstance(file, str):
            self.file = Path(file)

        if not self.file.is_file():
            raise FileNotFoundError(f"The file {self.file} does not exist.")

        self.run_type = run_type

        try:
            self.global_info = self.get_global_info(
                run_type=self.run_type, filename=self.file
            )
        except ValueError:
            raise ValueError(
                "-------------------------------------------\n"
                "Cannot parse the CP2K run_type information!\n"
                "Please check if you have provided an existing cp2k output file.\n"
                "If not, you can manually set the run_type parameter.\n"
                "Example:\n"
                "    Cp2kOutput(file='cp2k.out', run_type='md')\n"
            )
            sys.exit(1)

    def __str__(self):
        return f"Cp2kOutput(file={self.file}, run_type={self.run_type})"

    @staticmethod
    def get_global_info(run_type: str = "md", filename: Path | str = "cp2k.out"):
        """
        Get the global information from the CP2K output file.

        Parameters:
        -----------
        run_type: str
            The type of CP2K run. Options are 'md' for molecular dynamics and 'opt' for geometry optimization.
        filename: Path | str
            The path to the CP2K output file.

        Returns:
        --------
        dict
            A dictionary containing the global information.
        """
        if filename:
            global_info = parse_global_info(run_type=run_type.upper())
        elif run_type:
            global_info = GlobalInfo(run_type=run_type.upper())
        else:
            raise ValueError("Cp2k data does not know your run type!")
        return global_info
