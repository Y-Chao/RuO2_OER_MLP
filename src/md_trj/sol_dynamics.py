from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"


from typing import TYPE_CHECKING, Optional, Union

import numpy as np
from ase.atoms import Atoms
from ase.cell import Cell
from ase.geometry import get_distances
from ase.io import read
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.topology import tables
from scipy import constants

if TYPE_CHECKING:
    from MDAnalysis import Universe

"""
Statistical analysis of solvent dynamics at the solid-liquid interface.
"""


def modify_slab(slab: Union[str, Atoms], shift: float = 1e-4) -> Atoms:
    """
    Modify the slab to make the bottom of the slab close to the bottom of the box.
    """
    
    return new_slab

class SolMD(AnalysisBase):
    """
    Calulate the Solvent dynamic information of the system.
    """

    def __init__(
        self,
        universe: Union[str, Universe],
        slab_type: str = "asymmetry",
        topology: Optional[Union[str, Atoms]] = None,
        strict: bool = True,
        **kwargs,
    ) -> None:
        """
        Initialize the SolMD class.
        Parameters
        ----------
        universe : Universe
            The MDAnalysis Universe object.
        slab_type : str, optional
            The type of the slab, by default "asymmetry". Options are "asymmetry" and "symmetry".
        topology : str or Atoms,
            The topology file or Atoms object, by default None.
        strict : bool, optional
            Whether to use strict mode, by default True. If True, will raise error when no metal atoms found.
        kwargs : dict
            Additional arguments for defining the region, see `define_region` method.

        """
        if isinstance(universe, str):
            universe = Universe(universe)
        elif not isinstance(universe, Universe):
            raise TypeError("universe must be str or Universe object.")

        self.universe = universe
        self.slab_type = slab_type
        trajectory = universe.trajectory
        self.n_frames = len(trajectory)
        # Initialize the region
        self.topo, self.region = self.define_region(topology, **kwargs)
        super().__init__(trajectory, verbose=False)

        # Define the key species in the region
        self.O_ag = self.universe.select_atoms("name O")
        self.H_ag = self.universe.select_atoms("name H")
        self.M_ag = self.universe.atoms[self.choose_metal_atoms()]

        if self.universe.dimensions is not None:
            box = Cell.fromcellpar(self.universe.dimensions)
        else:
            box = self.topo.cell

        self.water_dic = {}

    def choose_metal_atoms(self):
        """
        Choose the metal atoms with atomic number > 18.
        Returns
        -------
        list
            The indices of the metal atoms in the universe.
        """
        elements = self.universe.atoms.elements
        element_numbers = [tables.SYMB2Z[element] for element in elements]
        metal_indices = [i for i, num in enumerate(element_numbers) if num > 18]
        if len(metal_indices) == 0:
            raise ValueError("No metal atoms found in the universe.")
        return metal_indices

    def define_region(self, topology: Optional[Union[str, Atoms]] = None, **kwargs):
        """
        Define the region of the simulation box.
        If the region is provied in the kwargs, use it directly.
        Otherwise, define the region based on the topology file or Atoms object.
        Parameters
        ----------
        topology : str or Atoms, optional
            The topology file or Atoms object, by default None.
        kwargs : dict
            Additional arguments for defining the region.
        Returns
        -------
        tuple
            A tuple containing the topology and the defined region.
        """

        # The region is defined by the topology or "region" or "cell" in kwargs
        if kwargs.get("region", None) is not None:
            region = del kwargs["region"]
        else: 
            if isinstance(topology, str):
                topo = read(topology)
            elif isinstance(topology, Atoms):
                topo = topology
            elif topology is None:
                symbols = self.universe.atoms.types
                positions = self.universe.select_atoms("all").positions
                if kwargs.get("cell", None) is None:
                    raise ValueError(
                        "Please provide the cell information when topology is None."
                    )
                topo = Atoms(
                    symbols=symbols,
                    positions=positions,
                    cell=Cell.fromcellpar(kwargs.get("cell", None)),
                )
            else:
                raise TypeError(
                    "topology must be str or Atoms object or provide the cell information."
                )

        # Define hegight of the region based on the positions of the metal atoms
        dim = kwargs.get("dim", 2)
        z_min = np.min(topo.positions[topo.numbers > 18][:, dim])
        z_max = np.max(topo.positions[topo.numbers > 18][:, dim])

        # Adjust the region to make the bottom of the slab close to the bottom of the box
        # There are three cases to consider:
        # 1. The slab is at the bottom of the box: Just minus the z_min
        # 2. The slab is in the middle of the box: Just minus the z_min
        # 3. The slab is separated at the top and bottom of the box: Adjust the positions of the upper bottom to the bottom

        if z_max > topo.cell.cellpar()[dim] / 2:
            tmp_z = z_max - topo.cell.cellpar()[dim]
            if tmp_z < z_min:
                z_max = z_min
                z_min = tmp_z
            else:
                raise ValueError(
                    "Please recheck the positions, to avoid the slab in the center of the box."
                )

        sol_region = [z_max, z_min + topo.cell.cellpar()[2]]
        topo.info["sol_region"] = sol_region

        return topo, sol_region

    def get_surface_metal_atoms(self, positions: np.ndarray):
        if self.slab_type == "asymmetry":
            z_max = np.max(self.M_ag.positions[:, 2])
            mask_z = self.M_ag[(z_max - 0.5)]
        elif self.slab_type == "symmetry":
            z_max = np.max(self.M_ag.positions[:, 2])
            z_min = np.min(self.M_ag.positions[:, 2])

            if z_max - z_min > self.topo.cell.cellpar()[2] / 2:
                tmp_z = z_max
                z_max = z_min
                z_min = tmp_z
            mask_z = self.M_ag[(z_max - 0.5) & (z_min + 0.5)]
        surface_positions = positions[self.M_ag.indices][mask_z]
        return surface_positions

    def identify_key_species(
        self,
        O_positions: np.ndarray,
        H_positions: np.ndarray,
        box: np.array | Cell,
        oh_cutoff: float = 1.3,
        **kwargs,
    ) -> np.ndarray:

        vec, dis = get_distances(p1=O_positions, p2=H_positions, cell=box, pbc=True)

        mask_dis = dis < oh_cutoff
        mask_CN = np.sum(mask_dis, axis=1) == 2

        water_dict = {}
        tolerance = kwargs.get("tolerance", 0.2)
        for i, is_water in enumerate(mask_CN):
            if is_water:
                pos_O = O_positions[i]
                pos_H = H_positions[mask_dis[i]]
                vec_OH = vec[i, mask_dis[i]]
                dis_OH = dis[i, mask_dis[i]]
                cos_HOH = np.degrees(
                    np.arccos(np.dot(vec_OH[0], vec_OH[1]) / (dis_OH[0] * dis_OH[1]))
                )

                # Double check the water molecules
                ## Check the angle between the two OH bonds
                if np.abs(cos_HOH - 104.5) > 15:
                    print(
                        f"Warning: Water molecule {i} is not valid. H-O-H angle: {cos_HOH:.2f}"
                    )
                    continue
                ## Check the positions of O
                # Should be noted that if the positions of O in periodic box, fix it later
                if (
                    pos_O[2] < self.region[0] + tolerance
                    or pos_O[2] > self.region[1] - tolerance
                ):
                    print(
                        f"Warning: Water molecule {i} is not valid. O position: {pos_O[2]:.2f}"
                    )
                    continue
                water_dict[i] = np.argwhere(mask_dis[i])
        return water_dict

    def calc_orientation(
        self,
        O_positions: np.ndarray,
        H_positions: np.ndarray,
        box: np.array | Cell,
        water_dict: dict,
    ) -> np.ndarray:
        """
        Calculate the orientation of the water molecules.
        """
        vec, dis = get_distances(p1=O_positions, p2=H_positions, cell=box, pbc=True)

        O_index = np.array([k for k in water_dict.keys()])
        H_index = np.array([v for v in water_dict.values()])

        rows = np.repeat(O_index, H_index.shape[1])
        cols = H_index.flatten()
        vec_OH = vec[rows, cols]
        dis_OH = dis[rows, cols]
        vec_OH = vec_OH.reshape((len(O_index), -1, 3))
        dis_OH = dis_OH.reshape((len(O_index), -1))

        return vec_OH

    def calc_position(
        self,
        O_positions: np.ndarray,
        H_positions: np.ndarray,
        water_dict: dict,
        water_ref: bool = False,
    ) -> np.ndarray:
        """
        Calculate the position of the water molecules.
        """
        O_index = np.array([k for k in water_dict.keys()])
        O_positions = O_positions[O_index]

        if water_ref:
            H_index = np.array([v for v in water_dict.values()])
            H_positions = H_positions[H_index]
            H_positions = H_positions.squeeze(axis=2)
            O_positions = O_positions[:, np.newaxis, :]

            all_positions = np.concatenate((O_positions, H_positions), axis=1)
            weights = np.ones_like(all_positions)
            O_weights = weights[:, 0] * 16
            H_weights = weights[:, 1:] * 1
            center_pos = np.average(all_positions, axis=1, weights=weights)
            return center_pos
        else:
            return O_positions

    def calc_possible_number(self) -> int:
        """
        Calculate the possible number of water molecules in the region.
        """
        z_min = self.region[0]
        z_max = self.region[1]
        num = np.sum(
            (self.O_ag.positions[:, 2] > z_min) & (self.O_ag.positions[:, 2] < z_max)
        )
        return num

    def _prepare(self):
        """
        Prepare the initialization of the analysis.
        """
        self.surface_area = np.linalg.norm(
            np.cross(self.topo.cell[0], self.topo.cell[1])
        )
        self.z1 = np.ones(self.n_frames) * self.region[0]
        self.z2 = np.ones(self.n_frames) * self.region[1]
        len_water = self.calc_possible_number()
        self.z_water = np.zeros((self.n_frames, len_water))
        self.orientation = np.zeros((self.n_frames, len_water, 3))
        self.angles = np.zeros((self.n_frames, len_water, 2, 3))
        self.water_coord = np.zeros((self.n_frames, len_water, 3))
        self.M_water = np.zeros(
            (self.n_frames, len_water, len(self.choose_metal_atoms()))
        )

    def _single_frame(self):
        ts_box = self.topo.cell

        positions = self._ts.positions
        O_positions = positions[self.O_ag.indices]
        H_positions = positions[self.H_ag.indices]
        # M_positions = positions[]

        new_water_dict = self.identify_key_species(
            O_positions, H_positions, ts_box, oh_cutoff=1.3
        )

        coord_water = self.calc_position(
            O_positions, H_positions, new_water_dict, water_ref=False
        )

        self.z_water[self._frame_index][: coord_water.shape[0]] = coord_water[:, 2]

        vec_OH = self.calc_orientation(O_positions, H_positions, ts_box, new_water_dict)

        _, dis_MO = get_distances(
            p1=coord_water,
            p2=positions[self.choose_metal_atoms()],
            cell=ts_box,
            pbc=True,
        )

        self.orientation[self._frame_index][: coord_water.shape[0]] = np.sum(
            vec_OH, axis=1
        )
        self.angles[self._frame_index][: coord_water.shape[0]] = vec_OH
        self.water_coord[self._frame_index][: coord_water.shape[0]] = coord_water
        self.M_water[self._frame_index][: coord_water.shape[0]] = dis_MO

    def _conclude(self):
        # self.z_water -= self.z1[:, np.newaxis]
        self.results.rho_water = self.calc_density_distribution(bin_width=0.1)
        self.results.orientation = self.calc_orientation_distribution(bin_width=0.1)
        self.results.ads_orientation = self.calc_KNN()

    def calc_KNN(
        self,
        M_OH_1: float = 2.5,
        M_OH_2: float = 4.5,
    ):
        ads_water_1st = []
        ads_water_2nd = []

        for i, dis_matrix in enumerate(self.M_water):
            min_RuO = np.min(dis_matrix, axis=1)
            O_Ru = (min_RuO <= M_OH_1) & (min_RuO > 0.5)
            lower_sec = self.water_coord[i][:, 2] < self.region[0] + 10
            mask_1st = O_Ru & lower_sec
            ads_water_1st.extend(
                SolMD.calc_theta(self.orientation[i][mask_1st], (0, 0, 1))
            )

        for i, dis_matrix in enumerate(self.M_water):
            min_RuO = np.min(dis_matrix, axis=1)
            O_Ru = (min_RuO > M_OH_1) & (min_RuO < M_OH_2)
            lower_sec = self.water_coord[i][:, 2] < self.region[0] + 10
            mask_2nd = O_Ru & lower_sec
            ads_water_2nd.extend(
                SolMD.calc_theta(self.orientation[i][mask_2nd], (0, 0, 1))
            )

        data_1 = self.calc_specific_orientation(np.array(ads_water_1st), bin_width=2)
        data_2 = self.calc_specific_orientation(np.array(ads_water_2nd), bin_width=2)
        return data_1, data_2

    def calc_density(self, n: float, volume: float, mol_mass: float = 18.015) -> float:
        rho = (n / constants.Avogadro * mol_mass) / (
            volume * (constants.angstrom / constants.centi) ** 3
        )
        return rho

    @staticmethod
    def calc_theta(vec: np.ndarray, ref_vec: tuple = (0, 0, 1)) -> float:
        """
        Calculate the angle between the vector and the reference vector.
        """
        return (
            np.arccos(
                np.dot(vec, ref_vec)
                / (np.linalg.norm(vec, axis=1) * np.linalg.norm(ref_vec))
            )
            / np.pi
            * 180
        )

    def calc_density_distribution(
        self,
        bin_width: float = 0.1,
    ) -> np.ndarray:
        """
        Calculate the density distribution of the water molecules.
        """
        z_min = self.region[0]
        z_max = self.region[1]

        bins = int((z_max - z_min) / bin_width)
        z_water = self.z_water.flatten()
        counts, bin_edges = np.histogram(z_water, bins=bins, range=(z_min, z_max))
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        n_water = counts / self.n_frames
        grid_volume = np.diff(bin_edges) * self.surface_area
        density = self.calc_density(n_water, grid_volume)
        return density, bin_centers

    def calc_orientation_distribution(
        self,
        bin_width: float = 0.1,
    ) -> np.ndarray:
        """
        Calculate the orientation distribution of the water molecules.
        """
        z_min = self.region[0]
        z_max = self.region[1]

        bins = int((z_max - z_min) / bin_width)

        hist, bin_edges = np.histogram(
            self.orientation.flatten(), bins=bins, range=(z_min, z_max)
        )
        hist = hist.astype(float)
        hist /= self.n_frames
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        return hist, bin_centers

    def calc_specific_orientation(
        self, ads_orientation: np.ndarray, bin_width: float = 2
    ) -> np.ndarray:
        """
        Calculate the specific orientation of the water molecules.
        """
        hist, bin_edge = np.histogram(
            ads_orientation, bins=90, range=(0, 180), density=True
        )
        bin_centers = bin_edge[:-1] + np.diff(bin_edge) / 2
        return np.array([bin_centers, hist]).T

    # def calc_mask_orientation(self,
    #                           bin_width: float = 0.1,
    #                           ) -> np.ndarray:
    #     z_min = self.region[0]
    #     z_max = self.region[1]
    #     bins = int((z_max - z_min) / bin_width)

    def calc_angle_distribution(
        self,
        bin_width: float = 2,
    ) -> np.ndarray:
        """
        Calculate the angle distribution of the water molecules.
        """
        angles = np.arccos(self.orientation.flatten()) / np.pi * 180 / self.n_frames
        angle_distribution, bin_edges = np.histogram(
            angles,
            bins=int(180 / bin_width),
            range=(0.0, 180.0),
            density=True,
        )
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        return angle_distribution, bin_centers


def main():
    example = "/Users/c_yang/OneDrive/nus/project/1_RuO2_stability/working/dataset/MLMD/RuO2_6L_sym_144_H2O/RuO2_6L_sym_144_H2O.xyz"
    universe = Universe(example)

    solmd = SolMD(universe)


if __name__ == "__main__":
    main()
