# JSBSim sounding rocket simulation

6-DOF flight simulation of a low-cost sounding rocket in [JSBSim](https://github.com/JSBSim-Team/jsbsim),
built from scratch for my master's thesis at the University of Queensland
(*Simulation Tools for Low-Cost Sounding Rockets*, 2026).

The rocket is a 2.5 m, 100 mm, ~21 kg single-stage vehicle on a Cesaroni M1670, aimed at 10,000 ft,
plus two-motor versions (cluster, series, two-stage). On top of the JSBSim model there are:

- a launch rail, air brake with apogee controller, dual-deploy recovery and a staging sequencer,
  all written as JSBSim flight-control channels
- a software-in-the-loop flight computer (sensor errors, 3- and 4-state EKF, apogee predictor)
- an independent 3-DOF point-mass code for verification
- Monte Carlo, one-at-a-time and Morris sensitivity studies
- validation against two flown rockets (Prometheus, Cavour)

## Setup

Python 3.11+.

```
pip install -r requirements.txt
```

No separate JSBSim install is needed; the `jsbsim` Python package is enough. Run everything from
the repo root. Don't run two scripts at once from the same folder: they share `data_output/`.

## Quick check

| Command | Should give |
|---|---|
| `python barrowman.py` | static margin 1.99 cal at launch, fin span 0.0463 m |
| `python sweep.py nominal --open` | apogee 10554.52 ft |
| `python sweep.py nominal` | apogee 10018.35 ft (brake on) |
| `python crosscheck.py verify` | 3-DOF apogee 10555 ft |
| `python stage_compare.py` | series 20549.7, cluster 20197.7, two-stage 20017.2 ft |
| `python sil.py nominal --ideal` | apogee 10013.9 ft |
| `cd kf && python test_kf.py` | ALL PASS |
| `cd val && python build_prometheus.py && python run_prometheus.py && python gps_check.py` | apogee 4186.5 m, +1.31 % vs GPS |
| `cd val && python build_cavour.py && python run_cavour.py && python cavour_check.py` | −1.64 % vs barometric record |

## Main studies

| What | Command | Time (1 core) |
|---|---|---|
| Launch angle × wind, 6-DOF vs point mass | `python anglesweep.py --root . --csv results_v6/anglesweep_v6.csv` | ~3 min |
| Recovery range of the brake | `python sweep.py recovery` | ~1 min |
| Brake size × gain map | `python sweep.py dmap` | ~3 min |
| Monte Carlo, 500 flights | `python sweep.py mcw 500 [--open \| --oracle]` | ~13 min |
| Flight computer in the loop | `python sil.py mc 500 --est kf4 --cd estimated` (`--ideal`, `--errors`, `--pad-temp-sd`, ...) | ~4 min |
| ±10 % sensitivity | `python sweep.py oat [--open]` | ~3 min |
| Morris screening | `python sweep.py morris 20 [--open]` | ~10 min |
| Tables from the result CSVs | `python batch2_tables.py --dir results_v6` | ~1 min |
| Brake loads, flap FE check | `python loads.py --root .`, `python flap_fea.py` | ~5 min |
| Staging | `python staging_seq.py coast\|wind\|tipoff\|mc --root . --tilt 1`, `python sweep_staging.py` | 1–30 min |

All Monte Carlo runs use seed 1, so they reproduce the numbers in `results_v6/`.
Use `-h` on any script for the options.

## Layout

```
aircraft/        JSBSim vehicle files (single stage, two-stage, cluster, series)
engine/          M1670 thrust table (from make_motor.py)
scripts/         launch.xml run script
kf/              Kalman filters used by sil.py
val/             Prometheus and Cavour validation models and flight data
atmos/           ERA5 temperature statistics
results_v6/      CSVs and figures behind the thesis
```

## Data

- Motor curves: Cesaroni RASP `.eng` files.
- Prometheus (Western Engineering Rocketry Team, 2022) and Cavour (PoliTo Rocket Team, 2023):
  vehicle data and flight logs from the data folder of the
  [RocketPy repository](https://github.com/RocketPy-Team/RocketPy), where the teams' permission is
  recorded. Only the data are used, not RocketPy itself.
- `atmos/atmos_stats.py` needs three ERA5 files from the same data folder
  (`EuroC_pressure_levels_reanalysis_2001-2021.nc`, `EuroC_single_level_reanalysis_2002_2021.nc`,
  `SpaceportAmerica_2018_ERA-5.nc`). They are not included here because of their size;
  the results are in `atmos/atmos_stats.json`. The same goes for the 2023 Spaceport America file used
  by `val/prometheus_inputs.py` (its output, `prometheus_inputs.json`, is included).

## Notes

- Units inside the JSBSim files are imperial, as JSBSim expects. `atmosphere/delta-T` is in Rankine.
- The air brake is modelled as drag only (no moment, no Mach dependence).
