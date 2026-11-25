#!/usr/bin/env python3
"""Recover the completed Q550 v18.53o server sendback without any ODE replay.

This is a post-processing-only recovery of the v18.53o package_sendback.py
PosixPath/string concatenation typo. Does not mutate the experiment root.
Python 3.9+ standard library only. All comments and messages in English.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import zipfile

ROOT_NAME = 'q550_full_10d_residual_valley_oriented_direction_boundary_source_transport_v18_53o_bundle'
EXPECTED_DECISIONS = {
    'AUTHENTIC_ORIENTED_10D_VALLEY_AND_BOUNDARY_GEOMETRY__NO_MATERIAL_HOLONOMY',
    'NUMERICAL_CLOSURE_SCREEN__PHYSICAL_RETURN_NOT_CERTIFIED',
}
REQUIRED_CSVS = {
    'fixed_width_oriented_transport.csv': 846,
    'local_valley_oriented_transport.csv': 846,
    'minimum_two_endpoint_orthogonality.csv': 6,
    'valley_curve_summary.csv': 12,
    'genuine_physical_ivps.csv': 6,
}
REQUIRED_FILES = [
    'results/report.json',
    'results/SCIENTIFIC_REPORT_SERVER.md',
    'results/derived_primitives.npz',
    'results/quadrature_pchip_refinement_sensitivity.csv',
]

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        while True:
            b = f.read(1024 * 1024)
            if not b: break
            h.update(b)
    return h.hexdigest()

def verify_manifest(root: Path, name: str, prefix: str) -> int:
    mpath = root / name
    if not mpath.is_file():
        raise RuntimeError(f'MISSING_{prefix}_MANIFEST {mpath}')
    obj = json.loads(mpath.read_text(encoding='utf-8'))
    manifest = obj['files']
    for rel, expected in manifest.items():
        rp = Path(rel)
        if rp.is_absolute() or '..' in rp.parts:
            raise RuntimeError(f'UNSAFE_{prefix}_MEMBER {rel}')
        p = root / rp
        if not p.is_file() or sha(p) != expected:
            raise RuntimeError(f'{prefix}_MANIFEST_MISMATCH {rel}')
    print(f'{prefix}_MANIFEST_PASS', len(manifest), flush=True)
    return len(manifest)

def count_csv_rows(path: Path) -> int:
    with path.open('r', encoding='utf-8', newline='') as f:
        return sum(1 for _ in csv.DictReader(f))

def verify_success(root: Path) -> dict:
    if root.name != ROOT_NAME:
        raise RuntimeError(f'WRONG_BUNDLE_ROOT {root.name} expected {ROOT_NAME}')
    for rel in REQUIRED_FILES:
        if not (root / rel).is_file():
            raise RuntimeError('MISSING_COMPLETED_RESULT ' + rel)
    report = json.loads((root / 'results/report.json').read_text(encoding='utf-8'))
    decision = report.get('scientific_decision')
    if decision not in EXPECTED_DECISIONS:
        raise RuntimeError('SCIENTIFIC_RUN_NOT_SUCCESSFUL ' + str(decision))
    if report.get('actual_Q550_ivps') != 10 or report.get('refined_Q550_ivps') != 4:
        raise RuntimeError('EXPECTED_10_PHYSICAL_IVPS_AND_4_REFINEMENTS')
    summary = report.get('summary', {})
    if summary.get('n_fixed_width_curve_rows') != 846 or summary.get('n_valley_follow_curve_rows') != 846:
        raise RuntimeError('INCOMPLETE_ORIENTED_CURVE_SUMMARY')
    for name, n_expected in REQUIRED_CSVS.items():
        p = root / 'results' / name
        if not p.is_file():
            raise RuntimeError('MISSING_CSV ' + name)
        actual = count_csv_rows(p)
        if actual != n_expected:
            raise RuntimeError(f'UNEXPECTED_CSV_ROWS {name} expected {n_expected} got {actual}')
    print('PHYSICAL_RESULTS_PRESENT 10 IVPS; 4 TIGHT; 846+846 CURVE ROWS', flush=True)
    return report

def verify_zip(zip_path: Path, root: Path, src_n: int, out_n: int) -> None:
    with zipfile.ZipFile(zip_path) as z:
        error = z.testzip()
        if error: raise RuntimeError('ZIP_CRC_FAIL '+error)
        names = z.namelist()
        expected_names = {str(Path(root.name) / p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
        if set(names) != expected_names:
            raise RuntimeError(f'ZIP_INCOMPLETE expected {len(expected_names)} entries got {len(names)}')
        for manifest_name, expected_count in [('SOURCE_MANIFEST.json',src_n), ('OUTPUT_MANIFEST.json',out_n)]:
            obj = json.loads(z.read(root.name+'/'+manifest_name))
            if len(obj['files']) != expected_count:
                raise RuntimeError('ZIP_MANIFEST_COUNT_CHANGED '+manifest_name)
            for rel, digest in obj['files'].items():
                if hashlib.sha256(z.read(root.name+'/'+rel)).hexdigest() != digest:
                    raise RuntimeError('ZIP_MANIFEST_HASH_MISMATCH '+rel)
        if json.loads(z.read(root.name+'/results/report.json'))['actual_Q550_ivps'] != 10:
            raise RuntimeError('ZIP_REPORT_IVPS_FAIL')
    print('ZIP_CRC_AND_BOTH_MANIFESTS_PASS',len(names),'members',flush=True)

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bundle-dir', type=Path, default=Path.cwd()/ROOT_NAME)
    ap.add_argument('--output', type=Path, default=None)
    ap.add_argument('--overwrite', action='store_true')
    args = ap.parse_args()
    root = args.bundle_dir.expanduser().resolve()
    if not root.is_dir():
        raise RuntimeError('BUNDLE_DIRECTORY_MISSING '+str(root))
    verify_success(root)
    n_source = verify_manifest(root,'SOURCE_MANIFEST.json','SOURCE')
    n_output = verify_manifest(root,'OUTPUT_MANIFEST.json','OUTPUT')
    out = (args.output.expanduser().resolve() if args.output else root.parent/(root.name+'_server_sendback.zip'))
    if out.parent != root.parent:
        print('NOTE: Output is outside the bundle root; this is allowed.',flush=True)
    if out.exists() and not args.overwrite:
        raise RuntimeError('OUTPUT_ALREADY_EXISTS '+str(out)+' (pass --overwrite to replace)')
    out.parent.mkdir(parents=True,exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=out.name+'.',suffix='.tmp',dir=str(out.parent))
    os.close(fd)
    temporary = Path(temporary)
    try:
        with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
            for p in sorted(root.rglob('*')):
                if p.is_file():
                    z.write(p, str(Path(root.name) / p.relative_to(root)))
        verify_zip(temporary,root,n_source,n_output)
        os.replace(temporary,out)
    finally:
        temporary.unlink(missing_ok=True)
    print('SEND_BACK_RECOVERED',out,flush=True)
    print('SHA256',sha(out),flush=True)
    print('NO_PHYSICAL_IVPS_EXECUTED',flush=True)
    return 0

if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print('RECOVERY_FAILED:',exc,file=sys.stderr)
        sys.exit(2)
