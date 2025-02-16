# MOORDS: Mooring Design and Simulation Package

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
- `moords/`: the core Python package for mooring design and simulation.
- `notebooks/`: contains Jupyter notebooks provides full design and simulation process. The notebook `design_example.ipynb` provides step-by-step guidance with detailed instructions. The nootebook `design_fieldwork.ipynb` provides an example of more complex moorings designed for fieldwork.
- `tests/`: includes Python test scripts to validate the functionality of the `moords` package.
- `output/`: contains outputs from the Jupyter notebooks, including printed tables, figures, and exported moorings for use in MATLAB.

<!-- ## Installation
To install the `moords` package, follow these steps:

```bash
git clone https://github.com/XXXX/moords.git
cd moords
pip install -r requirements.txt
``` -->
## Basic usage

To get started, we recommend following the step-by-step instructions provided in the `notebooks/design_example.ipynb` notebook. The general workflow is as follows:

1. Initialize a mooring object: create a new mooring object and specify the path to the mooring element database (e.g., `data/mooring_elements/mooring_elements_data.csv`).
2. Design the mooring: start from the bottom and add in-line elements in series. For clamp-on elements, define their height either along a specific element or along the entire mooring.
3. Inspect the design: use dataframes and design plots to review and validate the completed mooring design.
4. Generate reports: create an overview of the design, including serial numbers, and convert design figures into PDF format. Generate back-deck documents as PDFs for further documentation.
5. Export to MATLAB: export the mooring design for more detailed plotting and analysis in MATLAB.
6. Simulate the mooring: simulate the mooring under various flow conditions, either in a single instance or over a time series, and analyze the results.
