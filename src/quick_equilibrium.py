#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

"""
This module is to run equilibrium simulation quickly, to advoid the force field setting.
It is recommended to use uMLIP.
"""

import os

from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


def quick_equilibrium(
    atoms,
    model,
    temperature=300,
    pressure=1.0,
    nsteps=1000,
    timestep=1.0,
    ensemble="NVT",
    thermostat="Nose-Hoover",
    barostat="Nose-Hoover",
    device="cpu",
):
    """
    Run quick equilibrium simulation using uMLIP potential.

    1. NVT ensemble:
        - Langevin dynamics
        - Nose-Hoover dynamics
    2. NPT ensemble:
        - Berendsen

    Parameters
    ----------
    atoms : ASE Atoms object
        The input atomic structure.
    model : str
        The uMLIP model name or path.
    temperature : float, optional
        The temperature in Kelvin. Default is 300 K.
    pressure : float, optional
        The pressure in atm. Default is 1.0 atm.
    nsteps : int, optional
        The number of MD steps. Default is 1000.
    timestep : float, optional
        The time step in fs. Default is 1.0 fs.
    ensemble : str, optional
        The ensemble type: "NVT" or "NPT". Default is "NPT".
    thermostat : str, optional
        The thermostat type: "Nose-Hoover" or "Langevin". Default is "Nose-Hoover".
    barostat : str, optional
        The barostat type: "MelchionnaNPT" or "Berendsen". Default is "MelchionnaNPT".

    Returns
    -------
    atoms : ASE Atoms object
        The equilibrated atomic structure.
    """

    # if not os.getenv("ORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"):
    #     os.environ["ORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "0"

    if model == "uma":
        UMA_MLIP = os.getenv("UMA_MLIP")
        if UMA_MLIP == "":
            raise ValueError(
                "Please set the UMA_MLIP path by 'export UMA_MLIP=<your_path>'"
            )

        from fairchem.core import FAIRChemCalculator
        from fairchem.core.units.mlip_unit import load_predict_unit

        model = load_predict_unit(UMA_MLIP, device=device)
        calculator = FAIRChemCalculator(model, task_name="oc20", device=device)

        atoms.calc = calculator

    elif model == "mace":
        MACE_MLIP = os.getenv("MACE_MLIP")
        if MACE_MLIP == "":
            raise ValueError(
                "Please set the MACE_MLIP path by 'export MACE_MLIP=<your_path>'"
            )

        from mace.calculators import mace_mp

        calculator = mace_mp(model=MACE_MLIP, device=device)
        atoms.calc = calculator

    else:
        raise ValueError("Currently only support 'uma' and 'mace' models.")

    # Set the initial velocities corresponding to T=300K from Maxwell Boltzmann Distribution
    MaxwellBoltzmannDistribution(atoms, temperature_K=300)

    if ensemble == "NVT":
        if thermostat == "Langevin":
            from ase.md.langevin import Langevin

            # Langevin dynamics
            dyn = Langevin(
                atoms=atoms,
                timestep=timestep * units.fs,
                temperature_K=temperature,
                friction=0.02,
            )
            dyn.run(nsteps)
            return atoms

        elif thermostat == "Nose-Hoover":
            from ase.md.nose_hoover_chain import NoseHooverChainNVT

            # Nose-Hoover dynamics
            dyn = NoseHooverChainNVT(
                atoms=atoms,
                timestep=timestep * units.fs,
                temperature_K=temperature,
                tdamp=100 * timestep * units.fs,
                tchain=3,
            )
            dyn.run(nsteps)
            return atoms
        else:
            raise ValueError(
                "Currently only support 'Nose-Hoover' and 'Langevin' thermostats."
            )

    elif ensemble == "NPT":
        if barostat == "Berendsen":
            from ase.md.nptberendsen import NPTBerendsen

            dyn = NPTBerendsen(
                atoms=atoms,
                timestep=timestep * units.fs,
                temperature_K=temperature,
                pressure_au=atm2au(pressure),
                compressibility_au=1 / pa2au(1e11),  # 1/Pa
            )
            dyn.run(nsteps)
            return atoms
        else:
            raise ValueError("Currently only support 'Berendsen' barostat.")
    else:
        raise ValueError("Currently only support 'NVT' and 'NPT' ensembles.")


def pa2au(pressure_pa):
    """Convert pressure from Pa to atomic unit (au).
    1 Pa = 1 J/m^3 = 1.6e-19 eV/m^3 = 1.6e-29 eV/Angstrom^3
    """
    return pressure_pa * (units.eV / units.J) * (units.m**3)  # au


def atm2au(pressure_atm):
    """Convert pressure from atm to atomic unit (au).
    Here, use eV/Angstrom^3 as the atomic unit of pressure.
    """
    atm = 101325  # Pa

    return pa2au(pressure_atm * atm)
