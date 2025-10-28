#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"


import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np
from ase.io import read
from MDAnalysis import Universe

from md_trj.sol_dynamics import SolMD


def slice_trajectory(start, stop, step, n_chunks):
    frame_indices = list(range(start, stop, step))
    total = len(frame_indices)
    chunk_size = total // n_chunks
    extras = total % n_chunks

    splits = []
    s = 0
    for i in range(n_chunks):
        e = s + chunk_size + (1 if i < extras else 0)
        chunk_indices = frame_indices[s:e]
        if chunk_indices:
            splits.append(slice(chunk_indices[0], chunk_indices[-1] + step, step))
        s = e
    return splits


def run_single_solmd(
    split: slice, universe_path: str, topology_path: str = None, **kwargs
):
    u = Universe(universe_path)
    if topology_path is None:
        topo = read(universe_path, index=0)
    else:
        topo = read(topology_path)

    # Perform solvation analysis
    solmd = SolMD(universe=u, slab_type="assymmetry", topology=topo)
    solmd.run(split.start, split.stop, split.step)
    return solmd, solmd.results


def parallel_run_solmd(
    universe_path: str, start: int, stop: int, step: int, num_procs: int, **kwargs
):
    u = Universe(universe_path)
    splits = slice_trajectory(start, stop, step, num_procs)
    print(splits)

    run_partial = partial(run_single_solmd, universe_path=universe_path)
    with ProcessPoolExecutor(num_procs) as executor:
        results = list(executor.map(run_partial, splits))
    return results


def main():

    if len(sys.argv) < 6:
        print(
            "Usage: python parallel_solMD.py <universe_path> <num_procs> <start> <stop> <step>"
        )
        sys.exit(1)
    universe_path = sys.argv[1]
    # topology_path = sys.argv[2]
    num_procs = int(sys.argv[2])

    

    results = parallel_run_solmd(
        start=int(sys.argv[3]),
        stop=int(sys.argv[4]),
        step=int(sys.argv[5]),
        universe_path=universe_path,
        topology_path=None,
        num_procs=num_procs,
    )

    density = []
    ads_orientation = []
    ads_orientation_2nd = []
    for solmd, result in results:
        density.append([result["rho_water"][1], result["rho_water"][0]])
        ads_orientation.append(result["ads_orientation"][0])
        ads_orientation_2nd.append(result["ads_orientation"][1])

    density = np.array(density)
    ads_orientation = np.array(ads_orientation)
    ads_orientation_2nd = np.array(ads_orientation_2nd)

    density_out = np.mean(density, axis=0)
    ads_orientation_out = np.mean(ads_orientation, axis=0).T
    ads_orientation_2nd_out = np.mean(ads_orientation_2nd, axis=0).T
    np.save("universe_density.npy", density_out)
    np.save("universe_ads_orientation.npy", ads_orientation_out)
    np.save("universe_ads_orientation_2nd.npy", ads_orientation_2nd_out)


if __name__ == "__main__":
    main()
