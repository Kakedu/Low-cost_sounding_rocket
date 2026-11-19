"""Plot the four newly recorded nominal traces; no simulations are run here."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os

HERE = Path(__file__).resolve().parent
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE / "data",
                        help="Trace CSV directory (default: data beside this script)")
    parser.add_argument("--manifest", type=Path, default=HERE / "trajectory_run_manifest.json",
                        help="Recorded event/trace manifest")
    parser.add_argument("--output-dir", type=Path, default=HERE / "generated",
                        help="Figure and audit destination (default: generated beside this script)")
    args = parser.parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.manifest.expanduser().resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "SUCCESS"
    records = {r["case"]: r for r in manifest["cases"]}
    styles = [
        ("simultaneous", "Simultaneous", "#2166ac", "-"),
        ("series", "Series", "#b77712", "--"),
        ("staged_0s", "Staged: 0 s coast", "#24836a", "-."),
        ("staged_5s", "Staged: 5 s coast", "#b94545", ":"),
    ]
    traces = {}
    inputs = {}
    for name, *_ in styles:
        path = data_dir/f"{name}_trajectory.csv"
        assert sha(path) == records[name]["csv_sha256"]
        rows = list(csv.DictReader(path.open(encoding="utf-8")))
        traces[name] = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
        inputs[path.name] = sha(path)
    end_time = math.ceil((max(r["events"]["apogee"]["observed_time_s"] for r in records.values())+2)/5)*5
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 10, "axes.titlesize": 10.5, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "legend.fontsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .7, "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(7.4, 6.5), layout="constrained")
    grid = fig.add_gridspec(3, 1, height_ratios=[2.1, 2.0, 1.55])
    altitude = fig.add_subplot(grid[0])
    velocity = fig.add_subplot(grid[1], sharex=altitude)
    events_ax = fig.add_subplot(grid[2])
    plot_data = {}
    for name, label, color, linestyle in styles:
        trace = traces[name]
        selected = trace["t_s"] <= end_time
        x = trace["t_s"][selected]
        h = trace["h_agl_ft"][selected]/1000
        v = trace["v_vertical_fps"][selected]
        altitude.plot(x, h, color=color, ls=linestyle, lw=1.65, label=label)
        velocity.plot(x, v, color=color, ls=linestyle, lw=1.65)
        peak = records[name]["events"]["apogee"]
        altitude.plot(peak["observed_time_s"], peak["h_agl_ft"]/1000,
                      marker="o", ms=3, color=color, zorder=4)
        plot_data[name] = {"plotted_samples": int(selected.sum()),
                           "t_min_s": float(x.min()), "t_max_s": float(x.max()),
                           "h_min_kft": float(h.min()), "h_max_kft": float(h.max()),
                           "v_min_fps": float(v.min()), "v_max_fps": float(v.max())}
    altitude.set(title="(a) Height", ylabel="Height above ground (kft)", xlim=(0, end_time), ylim=(0, 23))
    altitude.legend(loc="upper left", ncols=2, frameon=False, handlelength=2.6, columnspacing=1.2)
    altitude.tick_params(labelbottom=False)
    velocity.set(title="(b) Vertical speed", ylabel="Vertical speed (ft/s)", xlabel="Time from ignition (s)")
    velocity.axhline(0, color="#777777", lw=.7)
    for ax in [altitude, velocity]:
        ax.grid(color="#dddddd", linewidth=.5)
        ax.set_axisbelow(True)
    plotted_events = {}
    for idx, (name, label, color, _) in enumerate(styles):
        lane = 3-idx
        record = records[name]
        ev = record["events"]
        positions = []
        for motor in [1, 2]:
            ignition = ev[f"motor{motor}_ignition_command"]["observed_time_s"]
            depletion = ev[f"motor{motor}_propellant_depleted"]["observed_time_s"]
            y = lane + (.105 if motor == 1 else -.105)
            events_ax.broken_barh([(ignition, depletion-ignition)], (y-.075, .15),
                                  facecolor=color if motor == 1 else "white", edgecolor=color, lw=.6,
                                  hatch=None if motor == 1 else "////", alpha=.85 if motor == 1 else 1.0)
            events_ax.plot(depletion, y, "o", color=color, mfc="white", mec=color, ms=4, zorder=4)
            if motor == 2:
                events_ax.plot(ignition, y, marker=">", color=color, mec="#222222", mew=.4, ms=5, zorder=5)
            positions.append({"motor": motor, "ignition_command_observed_s": ignition,
                              "propellant_depleted_observed_s": depletion})
        if "separation" in ev:
            # Offset the marker vertically so coincident 0s-coast separation
            # and second ignition remain separately visible at the same time.
            events_ax.plot(ev["separation"]["observed_time_s"], lane+.30, "D", color="#222222", ms=4, zorder=6)
        plotted_events[name] = positions
    events_ax.set(title="(c) Motor operation and separation", xlim=(-.15, 14), ylim=(-.55, 3.55),
                  xlabel="Time from ignition (s)", yticks=[3, 2, 1, 0],
                  yticklabels=[s[1] for s in styles], xticks=[0, 2, 4, 6, 8, 10, 12, 14])
    events_ax.tick_params(axis="y", length=0, labelsize=9)
    events_ax.grid(axis="x", color="#dddddd", linewidth=.5)
    event_legend = [
        Patch(facecolor="#888888", label="Motor 1"),
        Patch(facecolor="white", edgecolor="#888888", hatch="////", label="Motor 2"),
        Line2D([], [], marker=">", color="none", markerfacecolor="#888888", markeredgecolor="#222222", label="Motor 2 ignition", markersize=5),
        Line2D([], [], marker="o", color="none", markerfacecolor="white", markeredgecolor="#555555", label="Fuel depleted", markersize=4),
        Line2D([], [], marker="D", color="none", markerfacecolor="#222222", markeredgecolor="#222222", label="Separation", markersize=4),
    ]
    fig.legend(handles=event_legend, loc="outside lower center", ncols=3, frameon=False,
               fontsize=8.5, handlelength=1.4, columnspacing=1.3)
    png = output_dir/"Figure_7_1_staging_trajectories.png"
    pdf = output_dir/"Figure_7_1_staging_trajectories.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    caption = (
        "Nominal still-air trajectories for idealised simultaneous burn, series burn without mass release, "
        "and two-stage flight with 0 s or 5 s coast after separation. The launch tilt is 1 degree from vertical "
        "and the air brake is disabled. Panels (a) and (b) show height and vertical speed; panel (c) shows motor "
        "operation, second-motor ignition and separation. Motor-end symbols denote propellant depletion, "
        "not the thrust-below-1-lbf sequencing criterion. The staged cases apply a 0.2 s separation delay "
        "after the model's low-thrust detector, followed by the stated coast. Events are sampled at 1 ms. "
        "The configurations and apogees correspond to Tables 7.1 and 7.2; the unstaged references do not "
        "specify a buildable motor installation."
    )
    audit = dict(status="PASS_DATA_AND_RENDER_CREATED_PENDING_NATIVE_VIEW", figure=str(png),
                 figure_sha256=sha(png), pdf=str(pdf), pdf_sha256=sha(pdf), dpi=300,
                 figure_inches=list(fig.get_size_inches()), source_manifest_sha256=sha(manifest_path),
                 script_sha256=sha(__file__), input_csv_sha256=inputs,
                 axes_limits={"height": [list(altitude.get_xlim()), list(altitude.get_ylim())],
                              "vertical_speed": [list(velocity.get_xlim()), list(velocity.get_ylim())],
                              "events": [list(events_ax.get_xlim()), list(events_ax.get_ylim())]},
                 curve_data_summary=plot_data, event_markers=plotted_events,
                 max_abs_apogee_difference_ft=max(abs(r["apogee_difference_ft"]) for r in records.values()),
                 caption_suggestion=caption,
                 provenance_note="Four new deterministic reruns; this is not a recovered original trajectory record. "
                                 "Only plotting data windows are cropped; full near-ground CSVs are saved. "
                                 "No source model, accepted result, Monte Carlo draw, ignition gate or load factor was changed.")
    (output_dir/"figure_audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False), encoding="utf-8")
    (output_dir/"caption.txt").write_text(caption+"\n", encoding="utf-8")
    print(png)
    print("All four apogee differences:", [r["apogee_difference_ft"] for r in records.values()])
    plt.close(fig)


if __name__ == "__main__":
    main()
