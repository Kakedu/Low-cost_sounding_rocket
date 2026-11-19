"""Read-only numerical cross-check of the delivered ENGG7340 results.

Uses only Python's standard library. Reads existing CSV/JSON files, prints
the main results and exits nonzero on disagreement. Does not run JSBSim,
rewrite results, rebuild historical acceptance records or certify the model.
Put this script in the package root, or pass --root PATH.
"""

import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import sys


def load_json(path):
    with path.open(encoding="utf-8-sig") as stream:
        return json.load(stream)


def load_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def equal_number(actual, expected, message, tolerance=1e-9):
    require(math.isfinite(actual) and math.isfinite(expected), message + ": nonfinite number")
    require(math.isclose(actual, expected, rel_tol=1e-11, abs_tol=tolerance),
            f"{message}: found {actual:.12g}, expected {expected:.12g}")


def percentile(values, probability):
    """Linear interpolation, the accepted numpy.percentile(method='linear')."""
    values = sorted(values)
    require(bool(values), "Cannot take a percentile of an empty set")
    rank = (len(values) - 1) * probability / 100
    low = math.floor(rank)
    high = math.ceil(rank)
    return values[low] + (rank - low) * (values[high] - values[low])


def paired_result(first, second, count, allow_brake_difference=False):
    require([int(row["sample"]) for row in first[:count]] == list(range(count)),
            "First paired sample IDs are incomplete")
    require([int(row["sample"]) for row in second[:count]] == list(range(count)),
            "Second paired sample IDs are incomplete")
    inputs = ["isp", "cd0", "cn_alpha", "cp_x", "dry_mass", "drogue_cda",
              "main_cda", "launch_elev", "launch_az", "wind_fps", "wind_from", "dT_K"]
    if not allow_brake_difference:
        inputs.append("airbrake_dcd")
    for row_a, row_b in zip(first[:count], second[:count]):
        for name in inputs:
            equal_number(float(row_a[name]), float(row_b[name]),
                         f"Paired vehicle input {name}, sample {row_a['sample']}")
    a = [abs(float(row["apogee"]) - 10000) <= 100 for row in first[:count]]
    b = [abs(float(row["apogee"]) - 10000) <= 100 for row in second[:count]]
    first_only = sum(x and not y for x, y in zip(a, b))
    second_only = sum(y and not x for x, y in zip(a, b))
    discordant = first_only + second_only
    exact_p = (min(1.0, 2 * sum(math.comb(discordant, i)
                   for i in range(min(first_only, second_only) + 1)) / 2**discordant)
               if discordant else 1.0)
    return first_only, second_only, exact_p


def verify_control(root):
    summary = load_json(root / "rq2_revision2/final_control_summary.json")
    jobs = load_json(root / "study_protocol.json")["jobs"]
    primary = [job for job in jobs if job["kind"] in ("sil", "xml")]
    require(len(primary) == 63 and sum(job["n"] for job in primary) == 21900,
            "Formal control protocol is not 63 configurations / 21,900 outcomes")
    require(sum(job["n"] for job in primary if job["kind"] == "sil") == 16800,
            "SIL protocol count is not 16,800")
    require(sum(job["kind"] == "sil" for job in primary) == 48, "SIL set count is not 48")
    require(sum(job["n"] for job in primary if job["kind"] == "xml") == 5100,
            "True-state/open-loop protocol count is not 5,100")
    expected_names = {job["name"] for job in primary}
    actual_names = {path.name for path in (root / "results_revised").glob("mc_*.csv")}
    require(expected_names == actual_names, "Formal Monte Carlo CSV coverage differs from the protocol")
    data = {}
    for job in primary:
        name = job["name"]
        rows = load_csv(root / "results_revised" / name)
        rows.sort(key=lambda row: int(row["sample"]))
        require(len(rows) == job["n"], name + ": wrong number of flights")
        require([int(row["sample"]) for row in rows] == list(range(job["n"])),
                name + ": missing or duplicate sample IDs")
        heights = [float(row["apogee"]) for row in rows]
        require(all(math.isfinite(height) for height in heights), name + ": nonfinite apogee")
        target = float(job.get("target", 10000))
        hits = sum(abs(height - target) <= 100 for height in heights)
        accepted = summary["mc"][name]
        require(accepted["n"] == len(rows) and accepted["counts"]["hit"] == hits,
                name + ": accepted count disagrees with CSV")
        equal_number(accepted["rates"]["hit"], hits / len(rows), name + ": hit rate")
        if job["kind"] == "sil":
            require(all(0 < float(row["t_apogee"]) < float(row["t_land"]) < 599.9
                        for row in rows), name + ": invalid recorded flight endpoint times")
            require(all(math.isfinite(float(row[key])) for row in rows
                        for key in ("h_rms_m", "v_rms_ms")),
                    name + ": nonfinite full-ascent height/speed RMS")
        data[name] = rows

    cases = [
        ("Open loop", "mc_open_v6.csv", 47),
        ("True-state controller, design drag", "mc_closed_v6.csv", 369),
        ("True-state controller, known drag", "mc_oracle_v6.csv", 423),
        ("Three-state filter, ideal measurements", "mc_sil_kf_est_ideal_v6.csv", 397),
        ("Four-state filter, ideal measurements", "mc_sil_kf4_est_ideal_v6.csv", 397),
        ("Three-state filter, full sensor errors", "mc_sil_kf_est_v6.csv", 123),
        ("Four-state filter, full sensor errors", "mc_sil_kf4_est_v6.csv", 235),
        ("Four-state filter, known column temperature", "mc_sil_kf4_est_notemp_v6.csv", 385),
    ]
    print("Control: 63 formal configurations / 21,900 complete recorded outcomes")
    print("  48 SIL configurations / 16,800 + 15 open-loop/true-state configurations / 5,100")
    for label, name, expected_hits in cases:
        rows = data[name]
        hits = sum(abs(float(row["apogee"]) - 10000) <= 100 for row in rows)
        require(len(rows) == 500 and hits == expected_hits, label + ": differs from the thesis")
        print(f"  {label}: {hits}/500 = {100 * hits / 500:.1f}%")

    full = data["mc_sil_kf4_est_v6.csv"]
    height_rms = [float(row["h_rms_m"]) for row in full]
    require(all(math.isfinite(value) for value in height_rms), "Nonfinite main height RMS")
    for key, value in [("median", statistics.median(height_rms)),
                       ("p95", percentile(height_rms, 95)), ("max", max(height_rms))]:
        equal_number(value, summary["mc"]["mc_sil_kf4_est_v6.csv"]["h_rms_m"][key],
                     "Main full-ascent height RMS " + key)
    for key, expected in [("median", 15.0), ("p95", 43.0), ("max", 81.4)]:
        actual = {"median": statistics.median(height_rms),
                  "p95": percentile(height_rms, 95), "max": max(height_rms)}[key]
        equal_number(actual, expected, "Displayed thesis height RMS " + key, tolerance=0.05)
    print(f"  Four-state full-ascent height RMS: median {statistics.median(height_rms):.1f} m; "
          f"p95 {percentile(height_rms, 95):.1f} m; maximum {max(height_rms):.1f} m")
    global_h = max(float(row["h_rms_m"]) for name, rows in data.items()
                   if name.startswith("mc_sil_") for row in rows)
    global_v = max(float(row["v_rms_ms"]) for name, rows in data.items()
                   if name.startswith("mc_sil_") for row in rows)
    equal_number(global_h, summary["global_sil_tail_extrema"]["metrics"]["h_rms_m"]["value"],
                 "All-SIL height RMS maximum")
    equal_number(global_v, summary["global_sil_tail_extrema"]["metrics"]["v_rms_ms"]["value"],
                 "All-SIL speed RMS maximum")
    print(f"  Across all 48 SIL configurations: height RMS maximum {global_h:.1f} m; "
          f"speed RMS maximum {global_v:.1f} m/s. These tails remain in the statistics.")

    comparisons = [
        ("Exact temperature versus 3 K error", "mc_sil_kf4_est_notemp_v6.csv", "mc_sil_kf4_est_v6.csv", 500, 154, 4),
        ("Estimated versus design predictor drag", "mc_sil_kf4_est_v6.csv", "mc_sil_kf4_nom_v6.csv", 500, 46, 46),
        ("Known versus estimated predictor drag", "mc_sil_kf4_orc_v6.csv", "mc_sil_kf4_est_v6.csv", 500, 39, 14),
        ("Ideal fixed versus sampled brake, first 300", "mc_sil_kf4_est_ideal_fixbrake_v6.csv", "mc_sil_kf4_est_ideal_v6.csv", 300, 17, 0),
    ]
    for label, first, second, count, expected_first, expected_second in comparisons:
        only_a, only_b, p = paired_result(data[first], data[second], count,
                                        allow_brake_difference="_fixbrake_" in first)
        accepted = summary["paired"][first + " | " + second]
        require((only_a, only_b) == (expected_first, expected_second), label + ": discordant counts")
        require(accepted["n"] == count, label + ": accepted denominator")
        equal_number(p, accepted["mcnemar_exact_p"], label + ": exact paired p", tolerance=0.0)
        print(f"  {label}: n={count}, first-only/second-only={only_a}/{only_b}, exact p={p:.6g}")
    sampled = data["mc_sil_kf4_est_ideal_v6.csv"][:300]
    fixed = data["mc_sil_kf4_est_ideal_fixbrake_v6.csv"]
    require(sum(abs(float(row["apogee"]) - 10000) <= 100 for row in sampled) == 241,
            "Ideal sampled-brake prefix must be 241/300")
    require(sum(abs(float(row["apogee"]) - 10000) <= 100 for row in fixed) == 258,
            "Ideal fixed-brake prefix must be 258/300")
    print("  Ideal-brake comparison: 241/300 (80.3%) -> 258/300 (86.0%); not the 500-flight comparison")
    return summary


def verify_external_comparisons(root):
    base = root / "rq2_revision2/prometheus_temperature"
    comparison = load_json(base / "comparison_summary.json")
    inference = load_json(base / "temperature_inference.json")
    gps_reference = statistics.mean(item["gps_fit_peak_m"] for item in inference["flight_datasets"])
    cases = {case["name"]: case for case in comparison["cases"]}
    print("External flights (height definitions kept separate):")
    for name, expected_percent in [("regional_2022", 3.08), ("flight_geometric", 2.83)]:
        case = cases[name]
        rows = load_csv(base / "runs" / name / "model/data_output/prometheus.csv")
        apogee = max(float(row["h_agl_ft"]) for row in rows) * 0.3048
        error_percent = 100 * (apogee / gps_reference - 1)
        equal_number(apogee, case["apogee_geometric_m"], name + ": geometric apogee")
        equal_number(gps_reference, case["gps_reference_m"], name + ": fitted GPS reference")
        equal_number(error_percent, case["geometric_error_pct"], name + ": GPS percentage error")
        equal_number(error_percent, expected_percent, name + ": displayed thesis percentage", tolerance=0.005)
        print(f"  Prometheus {case['input_dT_K']:.2f} K: {apogee:.1f} m geometric AGL, "
              f"{error_percent:+.2f}% against fitted GPS {gps_reference:.1f} m "
              f"({case['role']})")
    check = load_json(root / "results_revised/validation/cavour_check.json")
    nominal = load_json(root / "results_revised/validation/cavour_dt_0.001/result.json")
    # Independent scalar form of the accepted pressure-altitude convention:
    # convert both geometric ASL heights to geopotential, then zero the
    # standard-atmosphere pressure-height difference at the 150 m pad.
    radius, pad = 6356766.0, 150.0
    geopotential = lambda h: radius * h / (radius + h)
    barometric = (geopotential(pad + nominal["apogee"]) - geopotential(pad)) * 288.15 / (288.15 + nominal["dT"])
    fitted_record = check["record_fit_m"]
    error_percent = 100 * (barometric / fitted_record - 1)
    equal_number(barometric, check["jsbsim_as_baro_m"], "Cavour model pressure height")
    equal_number(error_percent, check["E_pct"], "Cavour matched-barometric error")
    equal_number(error_percent, -2.03, "Cavour displayed thesis percentage", tolerance=0.005)
    print(f"  Cavour: {barometric:.1f} m model pressure height versus {fitted_record:.1f} m fitted "
          f"barometric record, {error_percent:+.2f}%")
    print("  GPS apex fits / the Cavour record fit are read from accepted observation JSONs; "
          "this helper does not refit the raw flight instruments.")


def verify_staging(root):
    result = root / "results_revised"
    nominal = load_csv(result / "staging_coast_tilt1_v6.csv")
    by_time = {float(row["coast"]): float(row["apogee"]) for row in nominal}
    gain = by_time[5.0] - by_time[0.0]
    equal_number(gain, 458.0, "Nominal still-air 1-degree five-second coast gain", tolerance=0.5)
    print("Staging (separate from the formal control denominator):")
    print(f"  Nominal 1-degree, still air: 0 -> 5 s coast gains {gain:.1f} ft")
    totals = 0
    for tilt, expected_lockouts in [(1, [0, 0, 19, 187]), (5, [0, 9, 150, 418])]:
        rows = load_csv(result / f"staging_mc_tilt{tilt}_v6.csv")
        require(len(rows) == 2000, f"Staging reference {tilt}: expected 2,000 trajectories")
        totals += len(rows)
        groups = {}
        for coast, lockouts in zip([0.0, 5.0, 10.0, 14.0], expected_lockouts):
            subset = sorted([row for row in rows if float(row["coast"]) == coast], key=lambda row: int(row["sample"]))
            require([int(row["sample"]) for row in subset] == list(range(500)),
                    f"Staging {tilt}-degree / {coast:g}s: incomplete sample IDs")
            require(sum(int(row["lockout"]) for row in subset) == lockouts,
                    f"Staging {tilt}-degree / {coast:g}s: lockout count")
            groups[coast] = subset
        print(f"  {tilt}-degree reference lockouts at 0/5/10/14 s: "
              + ", ".join(f"{n}/500 ({n/5:.1f}%)" for n in expected_lockouts))
        if tilt == 5:
            zero = statistics.mean(float(row["apogee"]) for row in groups[0.0])
            five = statistics.mean(float(row["apogee"]) for row in groups[5.0])
            ignited = [row for row in groups[5.0] if int(row["lit"]) == 1]
            equal_number(five - zero, -105.0, "Five-degree whole-sample mean change", tolerance=0.5)
            require(len(ignited) == 491, "Five-degree / five-second ignition subset is not 491")
            print(f"  5-degree / 5 s: whole-500 mean {zero:.1f} -> {five:.1f} ft, "
                  f"change {five-zero:+.1f} ft. Ignited-only mean "
                  f"{statistics.mean(float(row['apogee']) for row in ignited):.1f} ft uses n=491.")
    require(totals == 4000, "Staging trajectory count is not 4,000")
    print("  The four coast cases reuse each reference set's 500 vehicle inputs; "
          "lockouts remain in whole-sample means.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        verify_control(root)
        verify_external_comparisons(root)
        verify_staging(root)
    except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError, OverflowError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1
    print("PASS: stored primary results reproduce the checked thesis numbers.")
    print("The 4,500 force/attitude diagnostics, 7,000 historical-attribution flights, "
          "1,200 archived compensation trials and rejected repair attempts are not added to 21,900.")
    print("This is a read-only result check, not a new Monte Carlo run or a claim of hardware validation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
