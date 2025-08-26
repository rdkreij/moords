# MOORDS: Mooring Design and Simulation Package

> **Note**: The moords package is currently in the early stages of development and should be used with caution. While the core features are functional, there may still be bugs and incomplete features. User feedback is welcome as the package evolves.

The `moords` Python package provides advanced tools for the design, simulation, and visualization of mooring systems, made for oceanographic and marine engineering applications. It includes modules for designing mooring layouts, processing data, and simulation of the mooring in a flow. This package is specifically designed for fully submerged moorings.

The simulation algorithm is based on the Mooring Design & Dynamics MATLAB package by Richard K. Dewey ([MoorDyn](https://web.uvic.ca/~rdewey/mooring/moordyn.php)), which has been translated to Python. 

<!-- ## Features
- Mooring design: enables the configuration of both in-line and clamp-on mooring elements.
- Visualization: generates detailed tables and plots for mooring layouts, including serial numbers. The package also supports exporting mooring designs to MATLAB for additional detailed plotting.
- Automated reporting: converts design overviews into PDFs for streamlined documentation.
- Mooring simulation: simulates mooring behavior under various flow conditions, either as a single instance or over a time series.  -->

## Repository structure
- `data/`: contains the database with mooring element properties in the `mooring_elements` subfolder. It also includes simulated flow time series as NetCDF files in the `flow_series` subfolder, which are used as input for the mooring simulations.
- `matlab_plot/`: includes code for generating detailed plots of mooring designs. See `run_plot_code.m` for instructions.
- `src/`: the core Python package for mooring design and simulation.
- `notebooks/`: contains Jupyter notebooks provides full design and simulation process. The notebook `design_example.ipynb` provides step-by-step guidance with detailed instructions. The nootebook `design_fieldwork.ipynb` provides an example of more complex moorings designed for fieldwork.
- `tests/`: includes Python test scripts to validate the functionality of the `moords` package.
- `output/`: contains outputs from the Jupyter notebooks, including printed tables, figures, and exported moorings for use in MATLAB.

## Installation
The `moords` package can be installed through pip:

```bash
pip install git+https://github.com/rdkreij/moords
``` 

Alternative: set up a virtual environment using [Poetry](https://python-poetry.org/docs/) and run  
```bash
poetry install
```  


## Getting started
**Import required modules**
```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

import moords
```
**Define mooring element database:** create a dataframe of mooring elements with relevant properties.
```python
columns = ["name", "type", "buoyancy_kg", "length_m", "width_m", "diameter_m", "drag", "material", "comment"]
# Material index: 1: Steel, 2: Nylon, 3: Dacron, 4: Polyprop, 5: Polyethy, 6: Kevlar, 7: Aluminum, 8: Dyneema
rows = [
    ["3 Railway Wheels", "anchors", -1350.0, 0.55, 1.0, 0.0, 1.3, 1, ""],
    ["3/8 wire rope", "wires", -0.33, 1.0, 0.009, 0.0, 1.3, 1, ""],
    ["Zn Anode", "miscs", -0.05, 0.2, 0.025, 0.0, 1.3, 5, ""],
    ["41in ORE", "floats", 433.0, 1.04, 0.0, 1.04, 0.65, 1, ""],
    ["FLNTUSB", "miscs", -0.25, 0.31, 0.11, 0.0, 1.3, 5, ""],
    ["SBE39 T", "miscs", -0.25, 0.31, 0.11, 0.0, 1.3, 5, ""],
]
df_database = pd.DataFrame(rows, columns=columns).set_index("name")
```
**Build the mooring:** create a mooring system and add in-line and clamp-on elements.
```python
mooring = moords.Mooring(df_database, name="M1") 

# In-line elements
mooring.add_inline("3 Railway Wheels")
mooring.add_inline("3/8 wire rope", line_length=50, section="A")
mooring.add_clampon_along_inline(
    "Zn Anode", height_along_inline=4
)  # attach clamp-on to latest in-line: 3/8 wire rope
mooring.add_inline("41in ORE")

# Add clamp-on elements by height (can be done after adding all in-line elements)
mooring.add_clampon_by_height("FLNTUSB", serial="2927", height=20)
mooring.add_clampon_by_height("SBE39 T", serial="3921", height=30)
mooring.add_clampon_by_height("SBE39 T", serial="4263", height=35)

fig, ax = mooring.plot_design(label_rigging=True, line_ratio_plot=0.8)
plt.show()
```

**Simulate mooring behavior:** define water flow conditions and run the simulation.
```python
water_height = 60
z = np.linspace(0, water_height, 100)
ds_flow_instance = xr.Dataset(
    {
        "U": (["z"], z / water_height),
        "V": (["z"], np.zeros_like(z)),
        "W": (["z"], np.zeros_like(z)),
        "rho": (["z"], 1025 * np.ones_like(z)),
        "bottom_depth": water_height,
    },
    coords={"z": z},
)
instance = moords.SimulateInstance(mooring, ds_flow_instance)
instance.run()
fig, ax = instance.plot_instance()
plt.show()
```

## Basic usage

To get started, we recommend following the step-by-step instructions provided in the `notebooks/design_example.ipynb` notebook. The general workflow is as follows:

1. Construct a mooring element database containing the properties (e.g. length, buoyancy) of the inline and clampon elements. This can be constructed manually or loaded from a `.csv` (e.g., `data/mooring_elements/mooring_elements_data.csv`).
1. Initialize a mooring object: create a new mooring object and supply the mooring element database.
2. Design the mooring: start from the bottom and add in-line elements in series. For clamp-on elements, define their height either along a specific element or along the entire mooring.
3. Inspect the design: use dataframes and design plots to review and validate the completed mooring design.
4. Generate reports: create an overview of the design, including serial numbers, and convert design figures into PDF format. Generate back-deck documents as PDFs for further documentation.
5. Export to MATLAB: export the mooring design for more detailed plotting and analysis in MATLAB.
6. Simulate the mooring: simulate the mooring under various flow conditions, either in a single instance or over a time series, and analyze the results.
