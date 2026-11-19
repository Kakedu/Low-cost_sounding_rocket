"""New deterministic trajectory records for the supervisor's staging figure.

This does not alter the archived studies.  Each case uses byte-identical model
copies and the load/initialise/ignite/step sequence in stage_compare.fly.
The staged 5 s case applies the nominal staging_seq settings (1 deg reference
tilt, still air, k_stack=1, zero gate altitude, 20 deg tilt limit).
"""
from pathlib import Path
import argparse
import concurrent.futures
import csv
import datetime
import hashlib
import importlib
import json
import math
import os
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
def find_package_root():
    """Find the enclosing code package without relying on its folder name."""
    for candidate in [HERE, *HERE.parents]:
        if ((candidate / "stage_compare.py").is_file()
                and (candidate / "flight_checks.py").is_file()
                and all((candidate / name).is_dir() for name in ["aircraft", "engine", "scripts"])):
            return candidate
    raise FileNotFoundError("Code package root not found; pass --source-root explicitly")
DT = 0.001
CASES = {
    "simultaneous": dict(aircraft="myrocket_cluster", label="Simultaneous", coast_s=None,
                         ignite_both=True, engine_tanks=[0, 1],
                         reference="cluster, open loop"),
    "series": dict(aircraft="myrocket_series", label="Series", coast_s=None,
                   ignite_both=False, engine_tanks=[0, 1],
                   reference="series, open loop"),
    "staged_0s": dict(aircraft="myrocket2", label="Staged: 0 s coast", coast_s=0.0,
                      ignite_both=False, engine_tanks=[1, 0],
                      reference="two stage, open loop"),
    "staged_5s": dict(aircraft="myrocket2", label="Staged: 5 s coast", coast_s=5.0,
                      ignite_both=False, engine_tanks=[1, 0], reference=None),
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked_write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def worker(case_name, run_dir, snapshot):
    # Imports occur inside this fresh worker process.  Record and assert their
    # actual locations instead of relying solely on copied-file hashes.
    import jsbsim
    import numpy
    sys.path.insert(0, str(snapshot))
    sc = importlib.import_module("stage_compare")
    checks = importlib.import_module("flight_checks")
    assert Path(sc.__file__).resolve().is_relative_to(snapshot.resolve())
    assert Path(checks.__file__).resolve().is_relative_to(snapshot.resolve())
    case = CASES[case_name]
    start = time.perf_counter()
    f = jsbsim.FGFDMExec(str(run_dir))
    f.set_debug_level(0)
    checks.require(f.load_model(case["aircraft"]), "model load")
    f.set_dt(DT)
    checks.require(f.load_ic("reset", True), "initial conditions load")
    checks.require(f.run_ic(), "initial conditions execution")
    f.disable_output()
    f["fcs/target-apogee-ft"] = sc.NO_BRAKE_TARGET
    staged = case["coast_s"] is not None
    if staged:
        f["fcs/sep-delay-sec"] = 0.2
        f["fcs/coast-sec"] = case["coast_s"]
        f["fcs/gate-alt-ft"] = 0.0
        f["fcs/tilt-max-deg"] = 20.0
        f["fcs/use-inhibits"] = 1.0
        f["fcs/k-stack"] = 1.0
    initial = {k: f[k] for k in ["ic/theta-deg", "ic/psi-true-deg",
                                  "position/h-agl-ft", "inertia/weight-lbs",
                                  "atmosphere/wind-north-fps", "atmosphere/wind-east-fps"]}
    events = {"motor1_ignition_command": {"observed_time_s": 0.0}}
    if case["ignite_both"]:
        events["motor2_ignition_command"] = {"observed_time_s": 0.0}
    columns = ["t_s", "h_agl_ft", "v_vertical_fps", "mach", "qbar_psf", "weight_lb",
               "thrust0_lbf", "thrust1_lbf", "fuel0_lb", "fuel1_lb", "pitch_deg",
               "brake_position", "staged", "second_ignition_flag", "lockout"]

    def state():
        values = [f["simulation/sim-time-sec"], f["position/h-agl-ft"],
                  f["velocities/h-dot-fps"], f["velocities/mach"], f["aero/qbar-psf"],
                  f["inertia/weight-lbs"], f["propulsion/engine[0]/thrust-lbs"],
                  f["propulsion/engine[1]/thrust-lbs"], f["propulsion/tank[0]/contents-lbs"],
                  f["propulsion/tank[1]/contents-lbs"], f["attitude/theta-deg"],
                  f["fcs/airbrake-pos-norm"],
                  f["fcs/staged"] if staged else 0.0,
                  f["fcs/sust-lit"] if staged else
                  (f["fcs/second-burn"] if case_name == "series" else 1.0),
                  f["fcs/lockout"] if staged else 0.0]
        if not all(math.isfinite(v) for v in values):
            raise RuntimeError("Nonfinite trajectory state")
        return dict(zip(columns, values))

    def event(name, s, kind):
        if name not in events:
            events[name] = {"observed_time_s": s["t_s"], "h_agl_ft": s["h_agl_ft"],
                            "v_vertical_fps": s["v_vertical_fps"], "definition": kind}

    rows = [state()]
    peak = dict(rows[0])
    last_thrust_gt1 = [None, None]
    last_combined_thrust_gt1 = None
    lit = False
    step = 0
    while f["simulation/sim-time-sec"] < 200.0:
        if not lit and f["simulation/sim-time-sec"] >= sc.IGNITION_S:
            f["propulsion/engine[0]/set-running"] = 1
            f["fcs/throttle-cmd-norm"] = 1.0
            if case["ignite_both"]:
                f["propulsion/engine[1]/set-running"] = 1
                f["fcs/throttle-cmd-norm[1]"] = 1.0
            lit = True
        checks.run_checked(f)
        step += 1
        s = state()
        if s["h_agl_ft"] > peak["h_agl_ft"]:
            peak = dict(s)
        for engine, tank in enumerate(case["engine_tanks"]):
            if s[f"thrust{engine}_lbf"] > 0:
                event(f"motor{engine+1}_first_positive_thrust", s, "First observed thrust > 0 lbf")
            if s[f"thrust{engine}_lbf"] > 1:
                last_thrust_gt1[engine] = s["t_s"]
            if s[f"fuel{tank}_lb"] < 1e-10:
                event(f"motor{engine+1}_propellant_depleted", s, "First observed feeding-tank contents < 1e-10 lb")
        if s["thrust0_lbf"] + s["thrust1_lbf"] > 1.0:
            last_combined_thrust_gt1 = s["t_s"]
        if staged and s["staged"] > 0.5:
            event("separation", s, "First observed fcs/staged > 0.5; recorded at end of 1 ms step")
        if s["second_ignition_flag"] > 0.5:
            event("motor2_ignition_command", s, "First observed second-burn or sust-lit flag > 0.5")
        if staged and s["lockout"] > 0.5:
            raise RuntimeError("Unexpected nominal ignition lockout")
        # Store at 500 Hz, while event/maximum detection remains at every 1 ms step.
        if step % 2 == 0:
            rows.append(s)
        if s["t_s"] > 5.0 and s["h_agl_ft"] < 1.0:
            if not rows or rows[-1]["t_s"] != s["t_s"]:
                rows.append(s)
            break
    else:
        raise RuntimeError("Near-ground endpoint not reached by 200 s")
    events["apogee"] = {"observed_time_s": peak["t_s"], "h_agl_ft": peak["h_agl_ft"],
                        "v_vertical_fps": peak["v_vertical_fps"], "definition": "Maximum h-agl-ft over 1 ms integration states"}
    csv_path = run_dir.parent / f"{case_name}_trajectory.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, columns)
        writer.writeheader()
        writer.writerows(rows)
    reference_rows = list(csv.DictReader((snapshot/"reference_stage_compare_v6.csv").open()))
    if case["reference"] is not None:
        reference = next(float(r["apogee_ft"]) for r in reference_rows if r["configuration"] == case["reference"])
        reference_path = "results_revised/stage_compare_v6.csv"
    else:
        coast_rows = list(csv.DictReader((snapshot/"reference_staging_coast_tilt1_v6.csv").open()))
        reference = next(float(r["apogee"]) for r in coast_rows if float(r["coast"]) == 5.0)
        reference_path = "results_revised/staging_coast_tilt1_v6.csv"
    delta = peak["h_agl_ft"] - reference
    if abs(delta) > 1e-6:
        raise RuntimeError(f"Apogee differs from accepted nominal by {delta} ft")
    if max(abs(r["brake_position"]) for r in rows) > 1e-12:
        raise RuntimeError("Brake deployed in an unbraked architecture comparison")
    result = dict(case=case_name, configuration=case, status="SUCCESS", initial_conditions=initial,
                  integration_dt_s=DT, saved_trace_dt_s=0.002, full_near_ground_endpoint=rows[-1],
                  apogee_ft=peak["h_agl_ft"], reference_apogee_ft=reference,
                  apogee_difference_ft=delta, reference_result_path=reference_path,
                  events=events, last_individual_thrust_gt1_lbf_s=last_thrust_gt1,
                  last_combined_thrust_gt1_lbf_s=last_combined_thrust_gt1,
                  runtime_s=time.perf_counter()-start, csv_path=str(csv_path), csv_sha256=sha(csv_path),
                  actual_modules={"stage_compare": {"file": sc.__file__, "sha256": sha(sc.__file__)},
                                  "flight_checks": {"file": checks.__file__, "sha256": sha(checks.__file__)},
                                  "jsbsim": {"file": jsbsim.__file__, "version": jsbsim.__version__},
                                  "numpy": {"file": numpy.__file__, "version": numpy.__version__}},
                  actual_aircraft_path=str(run_dir/"aircraft"/case["aircraft"]/(case["aircraft"]+".xml")),
                  actual_aircraft_sha256=sha(run_dir/"aircraft"/case["aircraft"]/(case["aircraft"]+".xml")))
    checked_write_json(run_dir.parent/f"{case_name}_events.json", result)
    print(json.dumps({k: result[k] for k in ["case", "status", "apogee_ft", "apogee_difference_ft", "runtime_s"]}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path,
                        help="Code package root; default: nearest enclosing package directory")
    parser.add_argument("--output-dir", type=Path, default=HERE / "rerun_output",
                        help="New, nonexistent directory for an optional deterministic rerun")
    parser.add_argument("--workers", type=int, default=2,
                        help="Concurrent independent cases (default: 2)")
    parser.add_argument("--worker", choices=list(CASES))
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--snapshot", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.run_dir.resolve(), args.snapshot.resolve())
        return
    source = args.source_root.expanduser().resolve() if args.source_root else find_package_root()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise FileExistsError("--output-dir must not exist; preserve previous runs and recorded traces")
    output.mkdir(parents=True)
    snapshot = output/"source_snapshot"
    snapshot.mkdir()
    for folder in ["aircraft", "engine", "scripts"]:
        shutil.copytree(source/folder, snapshot/folder)
    for name in ["stage_compare.py", "staging_seq.py", "flight_checks.py", "motor_model.py",
                 "M1670.eng", "motor_build_info.json"]:
        shutil.copy2(source/name, snapshot/name)
    for name in ["stage_compare_v6.csv", "staging_coast_tilt1_v6.csv"]:
        shutil.copy2(source/"results_revised"/name, snapshot/("reference_"+name))
    bound = {q.relative_to(snapshot).as_posix(): sha(q) for q in snapshot.rglob("*") if q.is_file()}
    data = output/"data"
    data.mkdir(exist_ok=True)
    commands = []
    for name in CASES:
        run_dir = data/name
        run_dir.mkdir()
        for folder in ["aircraft", "engine", "scripts"]:
            shutil.copytree(snapshot/folder, run_dir/folder)
        (run_dir/"data_output").mkdir()
        commands.append([sys.executable, "-B", str(Path(__file__).resolve()), "--worker", name,
                         "--run-dir", str(run_dir), "--snapshot", str(snapshot)])
    start = time.perf_counter()

    def execute(command):
        name = command[command.index("--worker")+1]
        log = data/f"{name}_run.log"
        with log.open("w", encoding="utf-8") as stream:
            completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, cwd=output,
                                       env=dict(os.environ, JSBSIM_DEBUG="0"))
        if completed.returncode:
            raise RuntimeError(f"{name} failed; inspect {log}")
        print(f"Completed {name}", flush=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(execute, commands))
    results = [json.loads((data/f"{name}_events.json").read_text()) for name in CASES]
    assert all(sha(snapshot/name) == expected for name, expected in bound.items())
    manifest = dict(status="SUCCESS", recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    scope="Four NEW nominal deterministic reruns for a trajectory figure; not original trajectory logs and no new Monte Carlo samples.",
                    source_root=str(source), source_snapshot=str(snapshot), source_sha256=bound,
                    runner_path=str(Path(__file__).resolve()), runner_sha256=sha(__file__),
                    commands=commands, workers=args.workers, new_deterministic_cases=4, new_MC_flights=0,
                    runtime_s=time.perf_counter()-start, failures=0, cases=results,
                    event_time_convention="Events are first observed after each1ms run step; FCS command flags may refer to the preceding frame. No interpolation is presented as a more precise event measurement.",
                    endpoint="Full finite near-ground trajectory <1ft after5s before200s, matching stage_compare; plot uses ascent and early descent only.",
                    output_sha256={q.relative_to(output).as_posix(): sha(q) for q in data.glob("*") if q.is_file()})
    checked_write_json(output/"trajectory_run_manifest.json", manifest)
    print("ALL FOUR CASES PASS: nominal apogee difference <=1e-6ft", flush=True)


if __name__ == "__main__":
    main()
