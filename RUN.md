# Running the code

Run the commands from this folder. For new runs, make a working copy at a path
such as `D:/ENGG7340_work`. Some scripts replace model or result files.
Use an English path: the netCDF reader can fail on paths with Chinese characters.

## Check the saved files and results

```text
python verify_delivery.py
python verify_thesis_results.py
```

These commands only need Python's standard library. The first checks file
hashes and Python/XML syntax. The second recalculates the main thesis numbers
from the saved results. Neither runs a flight simulation.

## Install the packages

The simulations used Python 3.12.14. The package versions are in
`requirements-verified.txt`.

```text
python -m pip install -r requirements-verified.txt
```

JSBSim runs through its Python package; a separate executable is not needed.

## Model and filter checks

Run these in order, in a working copy:

```text
python barrowman.py
python verify_aero.py --json results_revised/aero_check.json
python verify_compiled_aero.py
python verify_reference_invariance.py
python verify_motor.py --json .reproduction/motor_check.json --scratch .reproduction/motor
python kf/test_kf.py
python crosscheck.py verify
python crosscheck_staging.py --jsb-csv results_revised/stage_compare_v6.csv
python check_rail.py
```

The compiled-force check reads the report from `verify_aero.py`.
Keep the `data_output/` directory, which is used by the JSBSim models.

For a single flight with ideal sensors:

```text
python sil.py nominal --est kf --cd estimated --ideal
python sil.py nominal --est kf4 --cd estimated --ideal
```

`kf` selects the three-state filter and `kf4` selects the four-state filter.
If `--est` is left out, the controller uses the true simulated states.

## Filter Monte Carlo runs

First prepare a new run folder without starting the flights:

```text
python rq2_revision2/rerun_sil_clean.py --prepare-only --destination D:/ENGG7340_sil_run
```

Then start the runs:

```text
python rq2_revision2/rerun_sil_clean.py --destination D:/ENGG7340_sil_run --workers 12 --resume
```

Use a new destination for a new study. This runs the 48 sensor-based
configurations, with 16,800 flight outcomes. It does not rerun the 15 open-loop
and true-state configurations. Their saved results are included here; the full
ZIP also has their original run records. Full studies can take hours.

## Physics, staging and flight comparisons

```text
python reproduce_physics.py --sandbox .reproduction/physics --workers 5
python val/prepare_prometheus_weather.py
python val/validation_mc.py
python atmos/atmos_stats.py
```

The first command runs the main physics and staging studies. Some supporting
checks were run separately; their commands and logs are in the full ZIP.
The validation command runs the flight comparisons and 200 uncertainty cases
for each vehicle.

## Draw figures from the saved results

```text
python batch2_tables.py --dir results_revised
python plot_nominal.py results_revised
python plot_results.py results_revised
python plot_margin.py --dir results_revised
python val/plot_prometheus_comparison.py
```

These commands write tables or figures from the existing data.

## Figures 3.1 and 7.1

These commands redraw the two figures from the included source and saved data:

```text
python figures/workflow/draw_workflow.py
python figures/staging/plot_staging_trajectories.py
```

Outputs go to each script's `generated/` folder. Both accept `--output-dir`.
Figure 7.1 uses the four nominal trajectories in `figures/staging/data/`; it
does not require another flight or Monte Carlo run. `figures/README.md`
explains the files and the optional deterministic rerun.
