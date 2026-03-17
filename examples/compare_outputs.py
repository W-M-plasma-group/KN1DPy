"""
Shared utility for comparing KN1DPy (Python) vs IDL outputs.
Called by the per-example check_outputs scripts.
"""

import numpy as np
from scipy.io import readsav


# Keys to skip (3D distribution functions with axis-ordering difference)
SKIP_KEYS = {'fh', 'fh2', 'fsh', 'eh_hist', 'si_hist'}

# Map Python npz key → IDL sav key (lowercase).
# Most are just lowercase; list exceptions here.
PY_TO_IDL = {
    'NuE': 'nue',
    'NuDis': 'nudis',
    'GaugeH2': 'gaugeh2',
    'PipeDia': 'pipedia',
    'WallH2': 'wallh2',
}


def _idl_key(py_key):
    """Return expected IDL key for a given Python key."""
    return PY_TO_IDL.get(py_key, py_key.lower())


def _rel_diff(py, idl):
    """Element-wise relative difference, guarded against divide-by-zero."""
    denom = np.maximum(np.abs(idl), 1e-30 * max(np.max(np.abs(idl)), 1.0))
    return np.abs(py - idl) / denom * 100.0


def _format_val(v):
    if np.isscalar(v) or v.ndim == 0:
        return f'{float(v):+.4e}'
    return f'{np.max(np.abs(v)):+.4e}'


def compare_file(py_npz, idl_sav_path, label):
    py   = np.load(py_npz)
    idl  = readsav(idl_sav_path)

    col = '{:<28s} {:>14s} {:>14s} {:>12s} {:>12s}'
    row = '{:<28s} {:>14s} {:>14s} {:>12.2f} {:>12.2f}'
    row_sc = '{:<28s} {:>14s} {:>14s} {:>12.2f} {:>12s}'
    sep = '-' * 84

    print()
    print('=' * 84)
    print(f'  {label}')
    print('=' * 84)
    print(col.format('Quantity', 'Python (|max|)', 'IDL (|max|)', 'Max %diff', 'Mean %diff'))
    print(sep)

    skipped = []
    missing = []

    for py_key in sorted(py.files):
        ik = _idl_key(py_key)
        if ik in SKIP_KEYS:
            skipped.append(py_key)
            continue
        if ik not in idl:
            missing.append(py_key)
            continue

        pv = py[py_key]
        iv = idl[ik]

        # Flatten / cast
        pv = np.asarray(pv, dtype=float).ravel()
        iv = np.asarray(iv, dtype=float).ravel()

        # If sizes differ, compare scalar Python value vs first IDL element
        if pv.size != iv.size:
            if pv.size == 1:
                iv = iv.ravel()[:1]
            else:
                missing.append(f'{py_key} (size mismatch: py={pv.size} idl={iv.size})')
                continue

        is_scalar = pv.size == 1

        if is_scalar:
            diff = _rel_diff(pv, iv)[0]
            print(row_sc.format(
                py_key,
                f'{float(pv[0]):+.4e}',
                f'{float(iv[0]):+.4e}',
                diff, '-'))
        else:
            diffs = _rel_diff(pv, iv)
            print(row.format(
                py_key,
                f'{np.max(np.abs(pv)):+.4e}',
                f'{np.max(np.abs(iv)):+.4e}',
                np.max(diffs),
                np.mean(diffs)))

    print(sep)
    if skipped:
        print(f'  Skipped (3-D arrays, axis ordering differs): {", ".join(skipped)}')
    if missing:
        print(f'  Not found in IDL output:')
        for m in missing:
            print(f'    {m}')


def run_comparison(case_name, py_dir, idl_dir, idl_prefix):
    print()
    print('#' * 84)
    print(f'#  KN1D output comparison: Python vs IDL   [{case_name}]')
    print('#' * 84)

    files = [
        ('KN1D_H',     f'{py_dir}/KN1D_H.npz',     f'{idl_dir}/{idl_prefix}.KN1D_H'),
        ('KN1D_H2',    f'{py_dir}/KN1D_H2.npz',    f'{idl_dir}/{idl_prefix}.KN1D_H2'),
        ('KN1D_input', f'{py_dir}/KN1D_input.npz',  f'{idl_dir}/{idl_prefix}.KN1D_input'),
        ('KN1D_mesh',  f'{py_dir}/KN1D_mesh.npz',   f'{idl_dir}/{idl_prefix}.KN1D_mesh'),
    ]

    for label, py_path, idl_path in files:
        compare_file(py_path, idl_path, label)

    print()
