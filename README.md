# Simulation tools for low-cost sounding rockets

Code, models and results for my master's thesis at the University of Queensland
(ENGG7340, *Simulation Tools for Low-Cost Sounding Rockets*, Haoran Du, October 2026).

The work builds a six-degree-of-freedom model of a low-cost sounding rocket in
[JSBSim](https://github.com/JSBSim-Team/jsbsim), driven from Python, and uses it to
answer three design questions:

1. **Model fidelity** – when is a point-mass model enough, and when are attitude dynamics needed?
2. **Apogee control** – how accurately can an air brake reach a 10,000 ft target when drag,
   motor performance and sensors are uncertain?
3. **Staging** – how do coast time, wind and ignition limits affect a two-stage design?

The reference vehicle is a 100 mm, single-stage rocket on a Cesaroni M1670 motor, with
two-motor variants (simultaneous, series and two-stage).

## Main results

| Question | Result |
| --- | --- |
| Fidelity | Point-mass and 6-DOF apogees agree within 0.075 % for near-vertical flight in still air, but differ by −28.65 % to +8.30 % in 20 mph wind. |
| Flight comparison | Prometheus (2022) geometric apogee +3.08 % against GPS with regional weather; Cavour (2023) −2.03 % against its barometric record. |
| Apogee control | Hit rate in the ±100 ft band rises from 9.4 % (no brake) to 73.8 % with true-state control. With all modelled sensor errors, the 4-state filter reaches 47.0 % and the 3-state filter 24.6 % (500 paired flights). |
| Staging | A 5 s coast raises apogee by 458 ft in still air, but at 5° launch tilt the full-sample mean falls once blocked ignitions are counted. |

The thesis gives the methods, assumptions and limits behind each number.

## What is in the repository

| Part | Files |
| --- | --- |
| JSBSim models | `aircraft/` (single stage, two-stage, simultaneous and series), `engine/`, `scripts/launch.xml` |
| Motor and aerodynamics | `motor_model.py`, `make_motor.py`, `barrowman.py`, `verify_*.py` |
| Flight computer (SIL) | `sil.py`, `kf/` (3- and 4-state extended Kalman filters, sensor and attitude models) |
| Parameter studies | `sweep.py` (Monte Carlo, one-at-a-time, Morris), `control_postprocess.py` |
| Point-mass cross-check | `crosscheck.py`, `crosscheck_staging.py` |
| Staging | `staging_seq.py`, `stage_compare.py`, `sweep_staging.py`, `staging_post.py` |
| Flight comparisons | `val/` (Prometheus and Cavour) |
| Brake loads | `loads.py`, `flap_fea.py`, `structure_screen.py` |
| Saved results | `results_revised/`, `rq2_revision2/` |
| Thesis figures | `thesis/figures/`, `figures/`, `THESIS_FIGURES.csv` |

`RUN.md` lists every command. `DATA_SOURCES.md` explains where the flight and weather
data come from.

## Setup

Python **3.12** is required (`sweep.py` uses syntax that Python 3.11 rejects).

```text
python -m pip install -r requirements-verified.txt
```

JSBSim runs through its Python package; no separate install is needed. Run the scripts
from the repository root and keep the `data_output/` folder.

Recalculate the main control results from the saved files (no simulation):

```text
python verify_thesis_results.py
```

This script stops at the Prometheus check, because the two Prometheus run traces
(`rq2_revision2/prometheus_temperature/runs/*/model/data_output/prometheus.csv`) were
not part of this copy. Rerun them with `python val/validation_mc.py` after adding the
flight logs listed under *Data not included*.

A single flight:

```text
python sweep.py nominal
python sil.py nominal --est kf4 --cd estimated --ideal
```

## Data not included

Some input files are not stored here because of their size or because they belong to
the RocketPy project. They are all in the `data` folder of the
[RocketPy repository](https://github.com/RocketPy-Team/RocketPy):

| File | Put it in | Needed for |
| --- | --- | --- |
| `2022-06-24-serial-5115-flight-0001-TeleMetrum.csv`, `2022-06-24-serial-6583-flight-0003-TeleMega.csv` | `val/data/` | Prometheus flight comparison |
| `EuroC_pressure_levels_reanalysis_2001-2021.nc`, `EuroC_single_level_reanalysis_2002_2021.nc`, `SpaceportAmerica_2018_ERA-5.nc` | `atmos/` | Weather statistics (`atmos/atmos_stats.py`); the results are already in `atmos/atmos_stats.json` |
| `spaceport_america_pressure_levels_2023_hourly.nc` | `val/data/` | Earlier 2023 proxy-weather case only |

Only the data are used; RocketPy itself is not run.

## Notes

- Units inside the JSBSim files are imperial, as JSBSim expects.
- Monte Carlo runs use fixed seeds, so they reproduce the saved results.
- Paired control cases reuse the same sampled vehicles: the 21,900 recorded flights
  are not 21,900 independent draws.
- `FILE_INDEX.csv` and `delivery_manifest.json` describe the full submission package,
  which also contained the large data files listed above.

## Licence

Apache 2.0 (see `LICENSE`). Flight data belong to the Western Engineering Rocketry Team
and the PoliTo Rocket Team, as published in the RocketPy repository.
