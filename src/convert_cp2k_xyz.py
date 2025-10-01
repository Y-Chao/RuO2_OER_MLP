#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

import argparse
import os
from glob import glob

import numpy as np
from ase import Atoms
from ase.io import read, write
from tqdm import tqdm

try:
    from cp2kdata import Cp2kOutput
except ImportError:
    raise ImportError("Please install cp2kdata package: pip install cp2kdata")

try:
    import dpdata
except ImportError:
    raise ImportError("Please install dpdata package: pip install dpdata")


def get_project_name(path: str) -> str:
    """
    Extract project name from a file path.

    Parameters:
    -----------
    path (str): File path.
    Returns:
    --------
    str: Project name.
    """
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if " GLOBAL| Project name" in line:
                project_name = line.split()[-1].strip()
                return project_name


def get_dump_xyz_freq(path: str) -> int:
    """
    Extract dump xyz frequency from a CP2K output file.

    Parameters:
    -----------
    path (str): Path to the CP2K output file.
    Returns:
    --------
    int: Dump xyz frequency.
    """
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "MD_PAR| Coordinates" in line:
                freq = int(line.split()[-2].strip())
                filename = str(line.split()[-1].strip())
                return filename, freq
    raise ValueError("Cannot find dump xyz frequency in the CP2K output file.")


def atoms_list(path: str, cell_info: np.ndarray, dump_freq: int = 1) -> list:
    """
    Extract atomic symbols from a CP2K output file.

    Parameters:
    -----------
    path (str): Path to the CP2K output file.
    Returns:
    --------
    list: List of atomic symbols.
    """
    atoms = read(path, index=":")
    if len(cell_info.shape) == 3:
        if int(cell_info.shape[0] / dump_freq) + 1 == len(atoms):
            for i, atom in enumerate(atoms):
                atom.set_cell(cell_info[i])
                atom.set_pbc([True, True, True])
        elif cell_info.shape[0] == 1:
            for i, atom in enumerate(atoms):
                atom.set_cell(cell_info[0])
                atom.set_pbc([True, True, True])
        else:
            raise ValueError(
                f"The number of cell info ({cell_info.shape[0]}) does not match the number of frames ({len(atoms)})."
            )
    elif len(cell_info.shape) == 2 and cell_info.shape == (3, 3):
        for i, atom in enumerate(atoms):
            atom.set_cell(cell_info)
            atom.set_pbc([True, True, True])
    else:
        raise ValueError("Cell info must be a (3, 3) or (n, 3, 3) array.")
    if isinstance(atoms, Atoms):
        return [atoms]
    return atoms


def SCF_convergence_list(path: str) -> bool:
    """
    Check if the SCF calculation in the CP2K output file has converged.

    Parameters:
    -----------
    path (str): Path to the CP2K output file.
    Returns:
    --------
    bool: True if SCF has converged, False otherwise.
    """
    converged = []
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if " SCF run " in line:
                if "converged" in line:
                    converged.append(True)
                else:
                    converged.append(False)
    return converged


def collect_cp2k_md(paths: list[str], output: str) -> list[Atoms]:
    """
    Collect MD trajectory from CP2K output files and save as XYZ format.

    Args:
        paths (list of str): List of paths to CP2K output files.
    """
    all_converged_atoms = []
    for path in paths:
        prefix_path = os.path.dirname(path)
        output_file = os.path.basename(path)
        cp2k_output = Cp2kOutput(path, prefix_path=prefix_path)
        num_frames = cp2k_output.num_frames
        scf_converge = SCF_convergence_list(path)
        if len(scf_converge) != num_frames:
            raise ValueError(
                f"[Warning] The number of SCF convergence records ({len(scf_converge)}) does not match the number of frames ({num_frames}) in {output_file}."
            )
        else:
            print(
                f"[Info] All SCF results equal to the number of frames in {output_file}."
            )

        atoms_converge = []
        energys = []
        forces = []
        stresses = []
        proj_name = get_project_name(path)
        traj_name, traj_freq = get_dump_xyz_freq(path)
        # print(traj_name, traj_freq)
        if proj_name is None:
            raise ValueError(
                f"[Error] Cannot find project name in {output_file}, please check your CP2K output file."
            )
        print(f"[Info] Project name: {proj_name}")
        atoms = atoms_list(
            os.path.join(prefix_path, traj_name), cp2k_output.all_cells, traj_freq
        )
        traj_frames = len(atoms)
        for t_frame in tqdm(
            range(traj_frames), total=traj_frames, desc="Processing frames"
        ):
            try:
                frame = t_frame * traj_freq
                if scf_converge[frame]:
                    if cp2k_output.has_force():
                        forces.append(cp2k_output.get_atomic_forces_list()[frame])
                    if cp2k_output.has_stress():
                        stresses.append(cp2k_output.get_stress_tensor_list()[frame])
                    energys.append(cp2k_output.get_energies_list()[frame])
                    tmp_atoms = atoms[t_frame].copy()
                    tmp_atoms.info["forces"] = forces[-1]
                    tmp_atoms.info["energy"] = energys[-1]
                    if cp2k_output.has_stress():
                        tmp_atoms.info["stress"] = stresses[-1]
                    atoms_converge.append(tmp_atoms)
                else:
                    print(
                        f"[Warning] SCF not converged in frame {frame} of {output_file}, skipping this frame."
                    )
            except:
                print(
                    f"[Error] Error processing frame {frame} of {output_file}, skipping this frame."
                )
                continue
        write(os.path.join(prefix_path, output), atoms_converge)
        all_converged_atoms.extend(atoms_converge)
    return all_converged_atoms


def parse_args():

    parser = argparse.ArgumentParser(
        description="Convert CP2K MD trajectory to XYZ format."
    )
    parser.add_argument(
        "-w",
        "--workdir",
        type=str,
        default=".",
        help="Working directory for output files.",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        type=str,
        default="log",
        help="Suffix of CP2K output files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="md_converged.xyz",
        help="Output XYZ file name.",
    )
    return parser.parse_args()


def main():
    print("[Usage] python convert_cp2k_xyz.py -h firstly to see the usage.")
    args = parse_args()

    paths = glob(os.path.join(args.workdir, f"*{args.suffix}"), recursive=True)
    atoms = collect_cp2k_md(paths, args.output)
    print(f"[Info] Collected {len(atoms)} converged frames from {len(paths)} files.")


if __name__ == "__main__":
    main()
