# ENGG7340 project files

This folder contains the main code, inputs and results for the thesis.
The full project is in `../03_Complete_Project/Haoran_Du_ENGG7340_Final_2026-10-29.zip`.
File paths follow those given in the thesis.

| Part | Files |
| --- | --- |
| Rocket models | `aircraft/`, `engine/` and `scripts/` contain the JSBSim models, motor inputs and launch settings. |
| Motor and aerodynamics | `motor_model.py`, `make_motor.py` and `barrowman.py` calculate the motor and aerodynamic inputs. The `verify_*.py` scripts check these calculations. |
| Flight validation | `val/` contains the Prometheus and Cavour scripts and flight data. Results are in `results_revised/validation/` and `rq2_revision2/prometheus_temperature/`. |
| Apogee control | `sil.py` runs the flight-computer cases. `kf/` contains the filters and sensor model. `sweep.py` runs the vehicle parameter studies. |
| Monte Carlo results | `results_revised/mc_*_v6.csv` contains the flight results. `rq2_revision2/final_control_summary.json` contains the control statistics. |
| Weather | `atmos/` contains the weather data and analysis. `atmosphere_model.py` supplies the atmosphere used by the simulation. |
| Staging | `staging_seq.py`, `stage_compare.py` and `sweep_staging.py` run the staging studies. `staging_post.py` calculates statistics from the saved flights. |
| Brake loads | `loads.py`, `flap_fea.py` and `structure_screen.py` calculate the loads and plate response. |

`THESIS_FIGURES.csv` links each thesis figure to its image, data and scripts.
The figure images are in `thesis/figures/`. `FILE_INDEX.csv` lists all the other
files in this folder.

`RUN.md` gives the package versions and commands. To check the saved thesis
results without rerunning flights, use `python verify_thesis_results.py`.

A few points matter when reading the results:

- Paired control cases reuse the same sampled inputs. The 21,900 recorded
  flight outcomes are not 21,900 independent input draws.
- Prometheus has a regional-weather prediction and a separate temperature
  calibration using its own flight measurements. `DATA_SOURCES.md` explains
  these comparisons and the assumed pressure-height conversion for Cavour.
- The staging results include flights where ignition was blocked. Keep those
  flights when calculating the mean apogee for the full sample.

The full ZIP contains the individual run folders and earlier studies omitted
here. `rq2_revision2/repair_source_v3/` is the code copy used by the filter
Monte Carlo launcher. Keep it unchanged if rerunning that study.

The current Word and PDF are in `thesis/`. Figure 3.1 source is in
`figures/workflow/`; Figure 7.1 source and its four nominal trajectories are
in `figures/staging/`. See `figures/README.md` for the plotting commands.

Appendix G prints the main model and analysis code with source line numbers.
`thesis/CODE_LISTINGS.csv` gives the file, printed range and thesis page for
each listing. Longer source files identified as extracts remain complete
in this folder. The thesis methods link to the matching appendix sections.

For the flap check, the thesis uses `sf_yield_beam` at the simulated load in
`results_revised/flap_fea_v6.json`. The saved file also includes scaled-load
and vibration results, which are outside the reported static comparison.

The plate geometry, material assumptions and mesh results are in Section F.4
and Table F.2. Sections 5.7 and 5.8 give the loads and the main static results.
