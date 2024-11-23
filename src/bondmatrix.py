#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'


# Fortran code may be faster for large matrix operations
# But it is hard to maintain and not easy to read.

class Bond:
    """
    Bond class is used to analysis and store the bond matrix information.
    """
    pass


def bondmatrix(
        const double[:, ::1] coords,
        const double[:, ::1] center_coord,
        const double cutoff,
        const np.int64_t[::1] pbc,
        const double[:, ::1] cell,
        const double tol=1e-8,
        const double min_r=1.0
    ):
    """This function is used to generate the bond matrix for the PBC structure.
    Parameters:
    -----------
    coords: np.ndarray
        The coordinates of all atoms in the system.
    center_coord: np.ndarray
        The coordinates of the center atom.
    cutoff: float
        The cutoff distance for the bond.
    pbc: np.ndarray
        The periodic boundary condition for the system.
    cell: np.ndarray
        The cell parameters for the system.
    tol: float
        The tolerance for the calculation.
    min_r: float
        The minimum distance for the bond.
    Returns:
    --------
    """

    if cutoff < min_r:
        raise ValueError('The cutoff distance is smaller than the minimum distance.')
    pass
    

    