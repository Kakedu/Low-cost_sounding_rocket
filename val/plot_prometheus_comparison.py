"""Plot regional and flight-informed Prometheus height comparisons."""
from pathlib import Path
import argparse
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent,
                        help="Project root; defaults to the parent of val/.")
    parser.add_argument("--output", type=Path,
                        help="Output PNG; defaults to results_revised/fig_prometheus_temperature_main.png.")
    args = parser.parse_args()
    root = args.root.resolve()
    data = root / "rq2_revision2/prometheus_temperature"
    output = args.output or root / "results_revised/fig_prometheus_temperature_main.png"
    comparison = json.loads((data / "comparison_summary.json").read_text(encoding="utf-8"))
    cases = {case["name"]: case for case in comparison["cases"]}
    regional = cases["regional_2022"]
    flight = cases["flight_geometric"]

    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for record, colour in [(regional, "#1764a0"), (flight, "#bd4d16")]:
        result = json.loads((data / "runs" / record["name"] / "result.json").read_text(encoding="utf-8"))
        pairs = [(t, h) for t, h in zip(result["t"], result["h"]) if t <= 40]
        ax.plot([t for t, h in pairs], [h for t, h in pairs], color=colour,
                label=f"JSBSim geometric: {record['input_dT_K']:.2f} K", lw=1.7)
    for name, filename, colour in [
        ("TeleMetrum", "2022-06-24-serial-5115-flight-0001-TeleMetrum.csv", "#222222"),
        ("TeleMega", "2022-06-24-serial-6583-flight-0003-TeleMega.csv", "#777777"),
    ]:
        with (root / "val/data" / filename).open(newline="", encoding="utf-8") as stream:
            rows = list(csv.reader(stream))
        header = [item.strip() for item in rows[0]]
        it = header.index("time")
        ih = header.index("height")
        ig = max(i for i, item in enumerate(header) if item == "altitude")
        pad = float(rows[1][ig])
        rows = [row for row in rows[1:] if len(row) == len(header) and 0 <= float(row[it]) <= 40]
        ax.plot([float(row[it]) for row in rows], [float(row[ig]) - pad for row in rows],
                color=colour, ls="--", lw=.9, label=name + " GPS")
        ax.plot([float(row[it]) for row in rows], [float(row[ih]) for row in rows],
                color=colour, ls=":", lw=.9, label=name + " barometric")
    ax.axvline(29.6, lw=.7, ls="--", color="#aaaaaa")
    ax.set(xlabel="Time from launch detection (s)", ylabel="Height above pad / initial GPS fix (m)",
           xlim=(0, 40), ylim=(0, 4500))
    ax.grid(alpha=.2)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(output.resolve())


if __name__ == "__main__":
    main()
