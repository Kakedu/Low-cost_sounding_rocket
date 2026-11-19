# Figure sources

The two plotting commands below read the supplied figure data. They do not run a flight simulation.

From the code package root, with NumPy and Matplotlib installed in the active Python environment:

```text
python figures/staging/plot_staging_trajectories.py
python figures/workflow/draw_workflow.py
```

Outputs go to each script's `generated` subdirectory. Both scripts also accept `--output-dir PATH`. Input defaults are relative to the script location, so the current working directory does not affect input discovery. The exact plotting environment used for the pixel comparison was Python 3.12.14, NumPy 2.5.4 and Matplotlib 3.11.2. Install the pinned plotting libraries if exact rendering is required:

```text
python -m pip install -r figures/requirements-figures.txt
```

## Figure 7.1

`figures/staging/data` contains the four full near-ground nominal trajectory CSVs and their exact event records. `trajectory_run_manifest.json` records the original deterministic runs and verifies the CSV hashes. The supplied paths in the historical records describe that earlier execution; the plotter does not use those absolute paths to locate inputs. The plots retain the original 1 ms event detection, 2 ms saved traces, apogee markers and separate propellant-depletion, second-ignition and separation markers. The plotted time window is shorter than the full supplied trajectories.

The four traces were recorded for the 28 October staging figure using the nominal parameters corresponding to Tables 7.1 and 7.2. They are deterministic cases, not additional Monte Carlo samples. `reference` contains the unchanged figure and original plotting audit.

A deterministic rerun is optional and is not needed to reproduce the figure. It additionally requires JSBSim 1.3.1 and the complete code package. The runner finds the nearest enclosing code package, or accepts `--source-root`. Its output directory must not already exist:

```text
python figures/staging/rerun_staging_trajectories.py --source-root . --output-dir .reproduction/staging_figure --workers 2
python figures/staging/plot_staging_trajectories.py --manifest .reproduction/staging_figure/trajectory_run_manifest.json --data-dir .reproduction/staging_figure/data --output-dir .reproduction/staging_figure/plots
```

The rerun copies aircraft, engine, launch inputs and the local driver modules into a new snapshot, retains the original numerical settings and rejects a case whose nominal apogee differs from the reference by more than 1e-6 ft. It never writes to the source model or supplied trace files.

## Figure 3.1

`figures/workflow/draw_workflow.py` generates `model_workflow.png` and its SVG from fixed drawing statements. It has no simulation-data inputs. Its `reference` image is unchanged.

`provenance/figure_sources_2026-10-29.json` records the source hashes, copied inputs and plotting checks.
