"""Recompute the frozen 48-set SIL study in a fresh, isolated project.

The accepted data, run histories, sources and package are never overwritten.
Without --resume, an existing destination is an error. --prepare-only creates
the isolated source/protocol layout without running the 16,800-flight study.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--destination', type=Path, default=ROOT / '.reproduction/revision2_sil')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--prepare-only', action='store_true')
    args = ap.parse_args()
    if args.workers < 1:
        raise ValueError('workers must be positive')
    destination = args.destination.resolve()
    if destination == ROOT.resolve() or destination == HERE.resolve():
        raise ValueError('Use a separate isolated destination')
    source = HERE / 'repair_source_v3'
    source_manifest = json.loads((source / 'source_manifest.json').read_text())
    expected = source_manifest['source_hashes']
    for relative, digest in expected.items():
        if sha(source / relative) != digest:
            raise ValueError('Frozen source changed: ' + relative)
    mini = destination / 'rq2_revision2'
    copied = {'source_manifest.json': source / 'source_manifest.json',
              'repair_protocol_v3.json': HERE / 'repair_protocol_v3.json',
              'run_accel_repair.py': HERE / 'run_accel_repair.py'}
    if args.resume:
        if not destination.is_dir():
            raise ValueError('No existing isolated run to resume')
        for relative, digest in expected.items():
            if sha(mini / 'repair_source_v3' / relative) != digest:
                raise ValueError('Isolated source mismatch: ' + relative)
        for name, original in copied.items():
            path = mini / 'repair_source_v3/source_manifest.json' if name == 'source_manifest.json' else mini / name
            if sha(path) != sha(original):
                raise ValueError('Isolated protocol/runner mismatch: ' + name)
    else:
        # mkdir(exist_ok=False) rejects existing data; this launcher never deletes.
        destination.mkdir(parents=True, exist_ok=False)
        mini.mkdir()
        shutil.copytree(source, mini / 'repair_source_v3')
        for name in ['repair_protocol_v3.json', 'run_accel_repair.py']:
            shutil.copy2(HERE / name, mini / name)
    command = [sys.executable, str(mini / 'run_accel_repair.py'), '--workers', str(args.workers), '--chunk-size', '100']
    if args.resume:
        command.append('--resume')
    record = dict(status='PREPARED', started_unix=time.time(), command=command,
                  launcher_sha256=sha(Path(__file__)), isolated_root=str(destination),
                  frozen_source_manifest_sha256=sha(source / 'source_manifest.json'),
                  source_hashes=expected,
                  protocol_sha256=sha(HERE / 'repair_protocol_v3.json'),
                  runner_sha256=sha(HERE / 'run_accel_repair.py'),
                  scope='Fresh48SIL/16800 study only; original 25 unaffected XML/nonMC jobs and the nominal trace are not recomputed by this command.')
    metadata = destination / ('launcher_resume_' + str(time.time_ns()) + '.json' if args.resume else 'launcher.json')
    metadata.write_text(json.dumps(record, indent=2))
    print('Isolated project:', destination, flush=True)
    print('Command:', subprocess.list2cmdline(command), flush=True)
    if not args.prepare_only:
        result = subprocess.run(command, cwd=destination, env=os.environ.copy())
        record.update(status='SUCCESS' if result.returncode == 0 else 'FAILED',
                      returncode=result.returncode, finished_unix=time.time())
        manifest = mini / 'repaired_runs_v3/repair_run_manifest.json'
        if manifest.is_file():
            record['execution_manifest_sha256'] = sha(manifest)
        metadata.write_text(json.dumps(record, indent=2))
        raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
