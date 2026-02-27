# Multiple Tidal Lagoons Hydrodynamics Model — UK

A 2D hydrodynamic model for simulating multiple idealised tidal range energy conversion (TREC) lagoons in the Severn Estuary and wider West UK coastal region, built on the [Thetis](https://thetisproject.org/) finite element ocean model.

---

## Scientific Background

The Severn Estuary has one of the largest tidal ranges in the world (up to ~14 m at spring tide), making it a prime candidate for tidal range energy extraction. Tidal lagoons — partially enclosed coastal impoundments fitted with low-head hydro turbines — can exploit this resource by generating electricity during both the ebb and flood phases of the tide.

While individual lagoon designs have been studied extensively, the **cumulative hydrodynamic impact** of deploying multiple lagoons simultaneously remains poorly understood. Modifying the tidal prism in one lagoon inevitably perturbs water levels and currents throughout the estuary, potentially affecting the performance of neighbouring lagoons and the broader coastal environment.

This model addresses that gap by simulating up to **seven idealised TREC lagoons** in the West UK simultaneously, allowing researchers to:

- Quantify inter-lagoon hydrodynamic interactions.
- Assess changes to tidal range, currents, and residual circulation.
- Evaluate aggregate power generation and capacity factors.
- Inform policy decisions on consenting multiple tidal range projects.

The work follows the consistent design approach described in Mackie *et al.* — *"Assessing impacts of tidal power lagoons of a consistent design"* — where each lagoon is parameterised identically to isolate location-dependent effects.

---

## Model Architecture

### Governing equations

The model solves the **2D depth-averaged shallow water equations** (SWE) using the [Thetis](https://thetisproject.org/) finite element solver, which is built on [Firedrake](https://www.firedrakeproject.org/). The time integration uses a Crank-Nicolson semi-implicit scheme (timestep Δt = 100 s) that is numerically stable for long tidal simulations without excessive diffusion.

Physical processes included:

| Process | Implementation |
|---|---|
| Bottom friction | Manning formulation (*n* = 0.03 s m⁻¹/³; spatially varying) |
| Tidal forcing | 8 constituents from TPXO (Q1, O1, P1, K1, N2, M2, S2, K2) |
| Coriolis | Constant-latitude *f*-plane (53 °N) |
| Wetting/drying | Enabled (α = 1.5, minimum depth = −10 m) |
| Viscosity sponge | Applied near open boundaries to absorb outgoing waves |

### Three-stage pipeline

```
0_preprocessing.py  →  1_ramp.py  →  2_run.py
```

| Stage | Duration | Purpose |
|---|---|---|
| **Preprocessing** | once | Build auxiliary fields: bathymetry, LAT, Manning friction, viscosity sponge |
| **Ramp-up** | 2 days | Spin up the model from rest; tidal forcing applied with a smooth ramp function |
| **Main simulation** | 30 days | Full operational run; lagoon turbines and sluices active |

### Lagoon operation state machine

Each lagoon is governed by a nine-mode state machine (`modules/lagoon_operation.py`):

- **Holding** — impoundment water level is allowed to build or drain.
- **Generating** — turbines operate when the head difference exceeds the minimum operating threshold (1 m).
- **Sluicing** — sluice gates open to equalise water levels.
- **Pumping** — turbines run in reverse to increase head before a generation phase.

Operational mode selection at each timestep depends on the instantaneous inner/outer water level difference (ΔZ) and the current operational phase (ebb or flood). A smooth ramp function prevents abrupt flow changes that could destabilise the finite element solver.

---

## Installation

### Recommended: Docker (simplest)

The project ships with a `Dockerfile` based on the official Firedrake image, which bundles Thetis, Firedrake, PETSc, and MPICH. No manual dependency installation is required.

```bash
# Build the image
docker build -t multiple-lagoons .

# Run the full simulation pipeline (outputs written to ./model_data)
docker run --rm -v "$(pwd)/model_data:/model_data" multiple-lagoons
```

### Manual installation

If you prefer to run outside Docker, you need:

1. **Firedrake** (includes PETSc, MPICH, UFL): follow the [Firedrake installation guide](https://www.firedrakeproject.org/download.html).
2. **Thetis**: install inside the Firedrake virtual environment:
   ```bash
   source firedrake/bin/activate
   pip install git+https://github.com/thetisproject/thetis.git
   ```
3. **Python dependencies** listed in `requirements.txt`:
   ```bash
   pip install -r requirements.txt
   ```
4. **UPTide** for tidal constituent predictions:
   ```bash
   pip install uptide
   ```

> **Note:** Firedrake must be activated before running any model scripts. All MPI-parallel execution uses `mpirun.mpich`.

---

## Input Data

Several large external datasets are required before running the model. These are **not included** in the repository due to size or licensing constraints, and must be obtained separately.

### Mesh files (provided in `inputs/`)

| File | Size | Description |
|---|---|---|
| `severn_outer_barrage.msh` | 5.6 MB | High-resolution Gmsh mesh of the Severn Estuary with lagoon boundaries |
| `ambient_UK_mesh10000.msh` | 2.0 MB | Coarser regional mesh covering the wider West UK shelf |

The mesh to use is selected in `inputs/simulation_parameters.py` via the `mesh_file` variable.

### Bathymetry (external — download separately)

| Dataset | Resolution | Usage |
|---|---|---|
| [DigiMap West UK](https://digimap.edina.ac.uk/) | 1 arc-second | Primary bathymetry source (high resolution) |
| [GEBCO](https://www.gebco.net/) | 15 arc-second | Fallback for areas outside DigiMap coverage |

Set the paths to these NetCDF files in `inputs/simulation_parameters.py`:

```python
bathymetry_paths = ["/path/to/digimap_west_uk.nc",
                    "/path/to/gebco_global.nc"]
```

### Tidal forcing (external — TPXO)

The model uses [TPXO](https://www.tpxo.net/) tidal constituent data. Request access from OSU and download:

- `gridES2008.nc` — TPXO grid file
- `hf.ES2008.nc` — tidal harmonic constants

Update the forcing paths in `inputs/simulation_parameters.py`:

```python
tpxo_grid = "/path/to/gridES2008.nc"
tpxo_data = "/path/to/hf.ES2008.nc"
```

Eight tidal constituents are extracted: **Q1, O1, P1, K1, N2, M2, S2, K2**.

### Gauge data (included)

`inputs/useful_gauges_BODC.csv` — 60+ British Oceanographic Data Centre (BODC) tidal gauge locations used as model validation points and in-model detectors.

---

## Configuration

All simulation settings are controlled from a single file: **`inputs/simulation_parameters.py`**.

### Key parameters

| Parameter | Default | Description |
|---|---|---|
| `mesh_file` | `west_uk_multiple_lagoons.msh` | Gmsh mesh to use |
| `dt` | `100` s | Model timestep |
| `ramp_time` | `2` days | Spin-up duration |
| `run_time` | `30` days | Main simulation duration |
| `n_manning` | `0.03` | Uniform Manning coefficient (overridden by spatially varying field) |
| `lat` | `53.0` °N | Latitude for Coriolis *f*-plane |
| `operation_mode` | `"two-way"` | Lagoon operation: `"ebb"`, `"ebb-pump"`, `"two-way"`, `"two-way-pump"` |
| `n_turbines` | `30` | Number of turbines per lagoon |
| `n_sluices` | `15` | Number of sluice gates per lagoon |
| `h_min` | `1.0` m | Minimum head for turbine generation |
| `rho` | `1025` kg m⁻³ | Seawater density |

### Lagoon identifiers

Seven lagoon sites are defined by their mesh boundary IDs:

| ID | Label | Location |
|---|---|---|
| 1 | SW | South Wales |
| 2 | CA | Cardiff |
| 3 | WA | Weston-super-Mare Area |
| 4 | CO | Colwyn Bay |
| 5 | LI | Liverpool |
| 6 | BL | Blackpool |
| 7 | SO | Solway Firth |

To run a subset of lagoons, comment out the unwanted entries in the `lagoon_ids` list.

### Turbine specification

```python
turbine_diameter   = 7.35   # m
turbine_capacity   = 20     # MW (per turbine)
turbine_efficiency = [0.95, 0.95]   # [generating, pumping]
sluice_area        = 100    # m² (total per lagoon)
sluice_Cd          = 1.0    # discharge coefficient
```

A template YAML alternative is provided in `inputs/template_parameters.yaml` for structured configuration management.

---

## Running the Model

### Full pipeline via shell script

```bash
bash sim.sh
```

`sim.sh` runs all three stages sequentially using 6 MPI processes:

```bash
mpirun.mpich -np 6 python 0_preprocessing.py
mpirun.mpich -np 6 python 1_ramp.py
mpirun.mpich -np 6 python 2_run.py
```

Adjust `-np 6` to match the number of available CPU cores.

### Running stages individually

```bash
# Stage 0: build auxiliary fields (run once per mesh/bathymetry change)
mpirun.mpich -np 6 python 0_preprocessing.py

# Stage 1: 2-day ramp-up — initialises ocean state from rest
mpirun.mpich -np 6 python 1_ramp.py

# Stage 2: 30-day production run
mpirun.mpich -np 6 python 2_run.py
```

> Stages must be run in order; each stage reads output from the previous one.

### Docker

```bash
docker run --rm \
  -v "$(pwd)/inputs:/inputs" \
  -v "$(pwd)/outputs:/outputs" \
  -v "/path/to/external/data:/data" \
  multiple-lagoons
```

### Visualisation

After the run completes, use the bundled plotting scripts:

```bash
# Single lagoon diagnostics (time series, power, head)
python plotting_single_lagoon.py

# Comparative plots across all lagoons
python plotting_multiple_lagoons.py
```

---

## Outputs

Outputs are written to `outputs/` and subdivided by pipeline stage.

### Preprocessed fields (`outputs/outputs_preproc/`)

| File | Description |
|---|---|
| `bathymetry2D.h5` | Interpolated bathymetry projected onto the mesh |
| `lat_field.h5` | Lowest Astronomical Tide correction field |
| `manning.h5` | Spatially varying Manning friction field |
| `viscosity_sponge.h5` | Near-boundary viscosity damping layer |

### Ramp-up outputs (`outputs/outputs_ramp/`)

HDF5 checkpoints of elevation and velocity fields at the end of the 2-day spin-up, used to initialise the main run.

### Main run outputs (`outputs/outputs_run/`)

**Field snapshots** (every 10 000 s ≈ 2.78 h):

| Variable | Description |
|---|---|
| `elev_2d` | Free-surface elevation (m) |
| `uv_2d` | Depth-averaged horizontal velocity (m s⁻¹) |

**Detector time series** (every timestep, 100 s):

- Tidal elevation at all BODC gauge locations and custom TRS detectors, stored as HDF5 datasets.

**Per-lagoon diagnostics** (every timestep):

Each lagoon `<ID>` produces `lagoon_<ID>.h5` containing:

| Variable | Unit | Description |
|---|---|---|
| `t` | s | Simulation time |
| `h_o` | m | Area-averaged outer (seaward) water level |
| `h_i` | m | Area-averaged inner (impoundment) water level |
| `DZ` | m | Head difference (h_o − h_i or h_i − h_o) |
| `P` | MW | Instantaneous power output |
| `E` | MWh | Cumulative energy generated |
| `Q_t` | m³ s⁻¹ | Turbine discharge |
| `Q_s` | m³ s⁻¹ | Sluice discharge |
| `mode` | — | Operational mode index (1–10) |

---

## Repository Structure

```
.
├── 0_preprocessing.py           # Stage 0: build bathymetry, friction, sponge fields
├── 1_ramp.py                    # Stage 1: 2-day spin-up simulation
├── 2_run.py                     # Stage 2: 30-day production simulation
├── plotting_single_lagoon.py    # Time series plots for a single lagoon
├── plotting_multiple_lagoons.py # Comparative plots across all lagoons
├── sim.sh                       # Shell script to run all three stages
├── Dockerfile                   # Container definition (Firedrake base image)
├── requirements.txt             # Additional Python dependencies
├── inputs/
│   ├── simulation_parameters.py   # Master configuration (edit this)
│   ├── template_parameters.yaml   # Alternative YAML config template
│   ├── useful_gauges_BODC.csv     # BODC tidal gauge coordinates
│   ├── extra_detectors_TRS.npy    # Custom TRS detector coordinates
│   ├── extra_detectors_TRS_names.npy
│   ├── severn_outer_barrage.msh   # Severn estuary mesh
│   ├── ambient_UK_mesh10000.msh   # Regional West UK mesh
│   ├── lat.txt                    # Latitude reference data
│   └── n_max_125.npy              # Manning coefficient field
├── modules/
│   ├── tools.py                   # LagoonCallback: per-timestep state updates
│   ├── lagoon_operation.py        # Nine-mode operational state machine
│   ├── input_barrages.py          # Lagoon initialisation and spec loading
│   ├── parameterisations.py       # Turbine hill charts, sluice discharge formulae
│   ├── input_0D.py                # 0D (box) model inputs
│   └── optimisation_functions.py  # Optimisation utilities
└── tools/
    ├── bathymetry.py              # Bathymetry interpolation to mesh
    ├── detectors.py               # BODC gauge and custom detector setup
    ├── thetis_support_scripts.py  # Coriolis, field initialisation, state export
    ├── tidal_forcing_ramp.py      # TPXO boundary condition updates
    ├── harmonic_analysis.py       # M2/S2 amplitude/phase decomposition
    ├── signal_processing.py       # Tidal signal reconstruction and range
    ├── gauge_analysis.py          # Gauge validation utilities
    ├── shear_stress.py            # Bed shear stress calculations
    ├── peaks.py                   # High/low water peak detection
    ├── tidal_amplitude.py         # Tidal amplitude calculations (parallel)
    ├── tidal_amplitude_serial.py  # Tidal amplitude calculations (serial)
    └── utm.py                     # UTM ↔ WGS84 coordinate conversion
```

---

## Citation

If you use this model in your research, please cite the reference publication:

> K Pappas, N.Q. Chien, I. Zilakos, L Beevers, A. Angeloudis, 2024. On
the economic feasibility of tidal range power plants. Proceedings of the
Royal Society A: Mathematical, Physical and Engineering Sciences. doi:
https://10.1098/rspa.2023.0867

> L. Mackie, S. C. Kramer, M. D. Piggott, and A. Angeloudis. Assessing impacts of tidal power lagoons of a consistent design. Ocean Engineering, 240:
109879, 2021. ISSN 0029-8018. doi: {https://doi.org/10.1016/j.oceaneng.2021.109879},
url: {https://www.sciencedirect.com/science/article/pii/S0029801821012282},


---

## References

- **TPXO tidal atlas**: Egbert, G.D. & Erofeeva, S.Y. (2002). Efficient inverse modeling of barotropic ocean tides. *Journal of Atmospheric and Oceanic Technology*, 19(2), 183–204.
- **Firedrake**: Rathgeber, F. *et al.* (2016). Firedrake: automating the finite element method by composing abstractions. *ACM Transactions on Mathematical Software*, 43(3), 24.
- **Thetis**: https://thetisproject.org/
- **BODC gauge data**: https://www.bodc.ac.uk/
- **GEBCO bathymetry**: https://www.gebco.net/

---

## Acknowledgements

This model was developed as part of PhD research at Cardiff University. The Thetis ocean model and Firedrake finite element framework are developed and maintained by their respective open-source communities. Bathymetry data were obtained from Ordnance Survey via Digimap and from the General Bathymetric Chart of the Oceans (GEBCO). Tidal constituent data are from the TPXO global tidal model (Oregon State University).

---

## Licence

This project is released for research use. Please contact the repository owner regarding any commercial applications.
