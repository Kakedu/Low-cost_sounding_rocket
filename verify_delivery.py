"""Verify the delivered package bytes using only the Python standard library."""
from pathlib import Path
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET

def file_sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def verify(root):
    root = root.resolve()
    manifest = json.loads((root / 'delivery_manifest.json').read_text(encoding='utf-8'))
    failures = []
    for name, expected in manifest['files_sha256'].items():
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f'Path outside package: {name}')
        actual = file_sha256(path) if path.is_file() else None
        if actual != expected:
            failures.append(dict(path=name, expected=expected, actual=actual))
    if failures:
        raise ValueError(json.dumps(failures, indent=2))
    py_count = xml_count = 0
    for name in manifest['files_sha256']:
        path = root / name
        if path.suffix == '.py':
            compile(path.read_bytes(), str(path), 'exec')
            py_count += 1
        elif path.suffix == '.xml':
            ET.parse(path)
            xml_count += 1
    return dict(status='PASS', verified_files=len(manifest['files_sha256']),
                python_sources_parsed=py_count, xml_sources_parsed=xml_count,
                scope='Delivered bytes and syntax only; flight and document checks are recorded separately.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    print(json.dumps(verify(args.root), indent=2))
