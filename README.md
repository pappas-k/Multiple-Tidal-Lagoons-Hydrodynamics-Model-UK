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
