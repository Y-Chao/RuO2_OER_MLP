# Extended Metal Oxide Interface (xMOI)
<!-- PROJECT LOGO -->
<div style="background-color: #D8BAE6; border-radius: 20px; padding: 0px;">
<p align="center">
<img src="docs/xMOI_logo.svg" alt="xMOI" style="width: 400px;">
</p>
</div>
Extended Metal Oxide Interface (xMOI) is a Python library designed to model, calculate and visualize the properties of metal (oxide) interfaces, particularly in the context of electrocatalysis.               
As our expectation, the library will provide tools for constructing surface/interface models, performing different calculations using various computational chemistry packages, and analyzing the results.

## Installation
You can install xMOI using pip:
```bash
git clone https://github.com/Y-Chao/xMOI.git
cd xMOI
pip install .
```
## Usage
### Build your models
1. Bulk structures
    * Set your `mp-api` key:
    ```bash
    export mp_api="your_api_key"
    ```
* Search your structures from Materials Project database:
    ```python
    from build_bulk import query_structure
    docs = query_structure(
            chemsystem=['Ru-O', 'Ru-O-H', 'Ru'],
            fields=['volume', 'material_id', 'symmetry']
            )
    ```
* Parse the structures and save them into a database:
    ```python
    from build_bulk import parse_bulk_data
    db = parse_bulk_data("bulk_example.db", docs)
    ```
* Prepare the DFT calculations:
    * Do calculations in place:
    ```python
    from build_bulk import prepare_calculation
    import os
    os.environ["VASP_PP_PATH"] = "/path/to/your/vasp/pseudopotentials/"
    prepare_calculation(db, xc='r2scan', calc_inplace=True)
    ```
    * Just prepare the input files:
    ```python
    from build_bulk import prepare_calculation
    import os
    os.environ["VASP_PP_PATH"] = "/path/to/your/vasp/pseudopotentials/"
    prepare_calculation(db, xc='r2scan', calc_inplace=False)
    ```

### Build your surface models
* Modify your configurations, which is placed in `examples/xMOI_configurations.toml`, for surface model construction, just comment the `solvation` line:
    ```toml
    [bulk_info]
    bulk_db = "bulk.db"      # [Mandatory] path to the ASE database file containing bulk structures
    sample = "RuO2"          # [Mandatory] name of the material sample
    crystal = "Rutile"       # [Mandatory] crystal structure
    xc = "pbe"               # [Mandatory] exchange-correlation functional, default: PBE

    [slab_info]
    miller_index = [0, 0, 1] # [Mandatory] miller index of the surface
    layers = 6               # [Mandatory] number of layers in the slab
    vacuum = 30.0            # [Mandatory] vacuum thickness in Angstrom
    layer_criteria = 0.1     # [Optional] minimum thickness of each layer in Angstrom
    terminations = []        # [Mandatory] List of possible terminations, or ["Ru", "O"]
    symmetry = true          # [Mandatory] whether to enforce symmetry when cleaving the slab
    fix = "bottom"           # [Optional] "center", "none"
    min_lattice = 10.0       # [Mandatory] minimum lattice parameter along x and y direction

    [interface_info]
    #solvation = "H2O"       # [Optional] solvation type, default: "H2O"
    num_sol = 0              # [Optional] number of solvent molecules, default: 0
    #sol_height = 10.0       # [Optional] height of the solvent box above the surface (in Angstrom), default: None
    surface_height = 1.0     # [Optional] add molecules above the surface by this height (in Angstrom)
    pH = 7                   # [Optional] pH value of the solution, default: 7
    ions = ["Cl-", "Na+"]    # [Optional] list of ions to add, e.g. ["Cl-", "Na+"]
    ions_number = [1, 1]     # [Optional] number of each ion, must match length of ions
    region = "bottom"        # [Optional] region to add solvent and ions, "top", "bottom", "both", default: "bottom"
    #seed = false            # [Optional] random seed for placing solvent and ions, default: None
    verbose = false          # [Optional] whether to print detailed information, default: false
    ```
* To run the surface generation code:
    ```bash
    python build_slab.py -i xMOI_configurations.toml
    ```
* The surface generation code is based on [pymatgen](https://pymatgen.org/) and [ASE](https://wiki.fysik.dtu.dk/ase/).

### Build your interface models
* Modify your configurations, which is placed in `examples/xMOI_interface_configurations.toml`, for interface model construction:
    ```toml
    [bulk_info]
    bulk_db = "bulk.db"          # [Mandatory] path to the ASE database file containing bulk structures
    sample = "RuO2"              # [Mandatory] name of the material sample
    crystal = "Rutile"           # [Mandatory] crystal structure
    xc = "pbe"                   # [Mandatory] exchange-correlation functional, default: PBE   

    [slab_info]
    miller_index = [0, 0, 1]     # [Mandatory] miller index of the surface
    layers = 6                   # [Mandatory] number of layers in the slab
    vacuum = 30.0                # [Mandatory] vacuum thickness in Angstrom
    layer_criteria = 0.1         # [Optional] minimum thickness of each layer in Angstrom
    terminations = ["Ru", "O"]   # [Mandatory] List of possible terminations, or ["Ru", "O"]
    symmetry = true              # [Mandatory] whether to enforce symmetry when cleaving the slab
    fix = "bottom"               # [Optional] "center", "none"
    min_lattice = 10.0           # [Mandatory] minimum lattice parameter along x and y direction

    [interface_info]
    solvation = "H2O"            # [Optional] solvation type, default: "H2O"
    num_sol = 20                 # [Optional] number of solvent molecules, default: 0
    sol_height = 10.0            # [Optional] height of the solvent box above the surface (in Angstrom), default: None
    surface_height = 1.0         # [Optional] add molecules above the surface by this height (in Angstrom)
    pH = 7                       # [Optional] pH value of the solution, default: 7
    ions = ["Cl-", "Na+"]        # [Optional] list of ions to add, e.g. ["Cl-", "Na+"]
    ions_number = [1, 1]         # [Optional] number of each ion, must match length of ions
    region = "bottom"            # [Optional] region to add solvent and ions, "top", "bottom", "middle", default: "bottom"
    #seed = false                # [Optional] random seed for placing solvent and ions, default: None
    verbose = false              # [Optional] whether to print detailed information, default: false
    ```
* Among the above parameters,
    * If you do not set the `num_sol`, it will be calculated based on the density of water at room temperature (1 g/cm3).
    * If you do not set the `sol_height`, it will be set to the vacuum thickness.
    * If you do not set the `surface_height`, it will be set to 1.0 Angstrom by default.
    * For `pH`, you have three options: 0, 7, and 14, which correspond to acidic, neutral, and basic conditions, respectively.
    * For `ions`, currently, only monovalent ions are supported.
    * For `region`, if you set it to "top", the solvent and ions will be added to the half of height of model; if you set it to "bottom", they will be added below the surface; if you set it to "middle", they will be added between above and below the surface.
* To run the interface generation code:
    ```bash
    python build_interface.py -i xMOI_interface_configurations.toml
    ```
* The interface generation code is based on [packmol](http://www.packmol.org/) and [ASE](https://wiki.fysik.dtu.dk/ase/).

## To Do
- [ ] To solve the reduction of water numbers when multiple ions added to the solvation box.
- [ ] To patch the surface models which utilize the `SlabGenerator` class in `pymatgen`.
- [ ] To patch the surface/interface models with adsorbates.
- [ ] Integrated with DFT packages, such as [VASP](https://www.vasp.at/), [CP2K](https://www.cp2k.org/), [GPAW](https://wiki.fysik.dtu.dk/gpaw/), etc. To provide the input files generation and output files parsing functions.
- [ ] Integrated with molecular dynamics packages, such as [LAMMPS](https://www.lammps.org/) etc. To provide the input files generation and output files parsing functions.
- [ ] Integrated with structure prediction packages, such as [AGOX](https://agox.gitlab.io/agox/index.html) etc. To provide the input files generation and output files parsing functions.
- [ ] Integrated with machine learning packages, such as [MACE](https://mace.readthedocs.io/en/latest/), [DeePMD-kit](https://deepmd.readthedocs.io/en/latest/), etc. To provide the input files generation and output files parsing functions.

## License
Distributed under the GPL-3.0 license. See `LICENSE` for more information.