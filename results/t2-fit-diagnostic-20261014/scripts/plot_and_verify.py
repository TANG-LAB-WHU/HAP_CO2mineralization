"""Verify the published T2 training evidence and optionally reproduce its plot.

Verification uses only the Python standard library. Plotting additionally needs
Matplotlib. No model, training dataset, runtime service or network is accessed.
The figures describe training fit; they do not establish independent validation.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import sys


SCHEMA = 'hap-T2-fit-diagnostic-publication/v1'
FIELDS = [
    'accepted_step', 'complete_objective_evaluations', 'elapsed_seconds',
    'force_component_RMSE_eV_A', 'max_abs_force_component_error_eV_A',
    'energy_RMSE_eV', 'energy_MAE_eV', 'max_abs_energy_error_eV',
    'energy_RMSE_meV_atom', 'energy_MAE_meV_atom', 'max_abs_energy_error_meV_atom',
    'objective', 'gradient_norm', 'gradient_max_abs',
]
THRESHOLDS = {
    'force_component_RMSE_eV_A_max': 0.01,
    'max_abs_force_component_error_eV_A_max': 0.05,
}
BRANCHES = [
    ('data/discarded-gpu-continuation.csv', 203, 2089, 2291, 'Discarded GPU continuation'),
    ('data/discarded-gpu-recovery.csv', 33, 2819, 2851, 'Discarded GPU recovery'),
    ('data/discarded-cpu1.csv', 142, 2995, 3136, 'Discarded CPU1'),
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_json(path):
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'Duplicate JSON key: ' + key)
            result[key] = value
        return result

    def reject_constant(value):
        raise ValueError('Nonfinite JSON value: ' + value)

    return json.loads(path.read_text(encoding='utf-8'),
                      object_pairs_hook=unique_pairs, parse_constant=reject_constant)


def bundle_file(root, relative):
    require(isinstance(relative, str) and relative and '\\' not in relative and ':' not in relative,
            'Artifact paths must be relative POSIX paths')
    parts = PurePosixPath(relative)
    require(not parts.is_absolute() and all(p not in ('.', '..') for p in parts.parts)
            and parts.as_posix() == relative, 'Unsafe artifact path: ' + relative)
    path = root.joinpath(*parts.parts)
    current = root
    for part in parts.parts:
        current = current / part
        require(not current.is_symlink(), 'Linked artifact is not allowed: ' + relative)
    require(path.is_file(), 'Missing artifact: ' + relative)
    require(root in path.resolve().parents, 'Artifact escaped bundle: ' + relative)
    return path


def file_digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def read_curve(path, start_step, end_step):
    rows = []
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames == FIELDS, 'Expected the original 14 CSV columns: ' + path.name)
        for line_number, raw in enumerate(reader, start=2):
            require(set(raw) == set(FIELDS) and all(raw[k] is not None for k in FIELDS),
                    'Incomplete or extra CSV fields at line ' + str(line_number))
            row = {}
            for key in FIELDS:
                try:
                    value = int(raw[key]) if key in FIELDS[:2] else float(raw[key])
                except (ValueError, TypeError):
                    raise ValueError('Invalid numeric value at line ' + str(line_number) + '/' + key)
                require(math.isfinite(value) and value >= 0,
                        'Nonfinite or negative curve value at line ' + str(line_number) + '/' + key)
                row[key] = value
            require(row['accepted_step'] == start_step + len(rows),
                    'Missing, duplicated or reordered accepted step: ' + path.name)
            if rows:
                require(row['complete_objective_evaluations'] > rows[-1]['complete_objective_evaluations']
                        and row['elapsed_seconds'] > rows[-1]['elapsed_seconds'],
                        'Curve evaluation/time counters must increase: ' + path.name)
            rows.append(row)
    require(len(rows) == end_step-start_step+1 and rows[-1]['accepted_step'] == end_step,
            'Curve row count or terminal step differs: ' + path.name)
    return rows


def verify_bundle(root):
    summary_path = bundle_file(root, 'summary.json')
    summary = load_json(summary_path)
    require(summary.get('schema') == SCHEMA, 'Unexpected publication summary schema')
    artifacts = summary.get('artifacts')
    require(isinstance(artifacts, dict) and artifacts, 'Missing artifact checksum manifest')
    verified_paths = {}
    for relative, item in artifacts.items():
        require(isinstance(item, dict) and type(item.get('bytes')) is int and item['bytes'] >= 0,
                'Invalid artifact byte count: ' + relative)
        digest = item.get('sha256')
        require(isinstance(digest, str) and len(digest) == 64
                and all(c in '0123456789abcdef' for c in digest), 'Invalid artifact SHA256: ' + relative)
        path = bundle_file(root, relative)
        require(path.stat().st_size == item['bytes'], 'Artifact byte count differs: ' + relative)
        require(file_digest(path) == digest, 'Artifact SHA256 differs: ' + relative)
        verified_paths[relative] = path.resolve()

    curve = summary['curve']
    require(curve.get('file') == 'data/curve.csv' and type(curve.get('rows')) is int
            and type(curve.get('accepted_updates')) is int
            and curve['rows'] == 10979 and curve['accepted_updates'] == 10978,
            'Expected the published 0..10978 retained curve')
    require(curve['file'] in verified_paths, 'Main curve is absent from the checksum manifest')
    rows = read_curve(verified_paths[curve['file']], 0, curve['accepted_updates'])
    require(len(rows) == curve['rows'], 'Summary/main-curve row count differs')
    require(summary.get('thresholds') == THRESHOLDS, 'Original force thresholds changed')
    metrics = summary['metrics']
    require(isinstance(metrics, dict) and set(FIELDS[3:5]).issubset(metrics),
            'Terminal force metrics are missing')
    for key, value in metrics.items():
        require(key in FIELDS[3:] and type(value) in (int, float) and math.isfinite(value),
                'Unknown or nonfinite terminal metric: ' + key)
        require(math.isclose(rows[-1][key], value, rel_tol=1e-12, abs_tol=1e-14),
                'Terminal metric differs from the final curve row: ' + key)
    first_pass = next((row['accepted_step'] for row in rows
                       if row[FIELDS[3]] <= THRESHOLDS['force_component_RMSE_eV_A_max']
                       and row[FIELDS[4]] <= THRESHOLDS['max_abs_force_component_error_eV_A_max']), None)
    require(first_pass == curve['accepted_updates'],
            'The final accepted step must be the first retained step meeting both force thresholds')

    declarations = summary['discarded_branches']
    require(isinstance(declarations, list) and len(declarations) == len(BRANCHES),
            'Expected exactly three discarded branches')
    require(all(isinstance(item, dict) for item in declarations), 'Invalid discarded-branch declaration')
    by_file = {item['file']: item for item in declarations}
    require(len(by_file) == len(BRANCHES) and set(by_file) == {b[0] for b in BRANCHES},
            'Discarded branch filenames are missing or duplicated')
    branches = []
    for relative, count, start, end, label in BRANCHES:
        item = by_file[relative]
        require(all(type(item.get(k)) is int for k in ('accepted_updates', 'start_step', 'end_step'))
                and (item['accepted_updates'], item['start_step'], item['end_step']) == (count, start, end),
                'Discarded branch count/range differs: ' + relative)
        require(relative in verified_paths, 'Discarded branch is absent from the checksum manifest')
        branch_rows = read_curve(verified_paths[relative], start, end)
        require(len(branch_rows) == count, 'Discarded branch row count differs: ' + relative)
        branches.append((label, branch_rows))
    report = dict(status='verified', schema=SCHEMA, artifacts_checked=len(verified_paths),
                  curve_rows=len(rows), accepted_updates=rows[-1]['accepted_step'],
                  first_passing_accepted_step=first_pass, discarded_branch_updates=[len(b[1]) for b in branches],
                  independent_validation=False, model_evaluations=0, training_only=True)
    protected_paths = set(verified_paths.values()) | {summary_path.resolve()}
    return rows, branches, report, protected_paths


def plot_curves(rows, branches, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42})
    figure, axes = plt.subplots(2, 2, figsize=(12, 8))
    panels = [
        (FIELDS[3], 'Force component RMSE', 'eV / Angstrom', 0.01),
        (FIELDS[4], 'Maximum absolute force component error', 'eV / Angstrom', 0.05),
        ('energy_RMSE_meV_atom', 'Energy RMSE (report only)', 'meV / atom', None),
        ('objective', 'Training objective', 'Weighted objective', None),
    ]
    series = [('Retained accepted sequence', rows, '#17694d', '-')]
    series += [(label, values, color, '--') for (label, values), color in
               zip(branches, ('#c77b20', '#7752a0', '#bd414a'))]
    for axis, (key, title, unit, threshold) in zip(axes.flat, panels):
        for label, values, color, style in series:
            axis.plot([row['accepted_step'] for row in values], [row[key] for row in values],
                      color=color, linestyle=style, linewidth=1.2, label=label)
        if threshold is not None:
            axis.axhline(threshold, color='#555555', linestyle=':', linewidth=1,
                         label='Training threshold ' + str(threshold))
        if key == 'objective' and all(row[key] > 0 for _, values, _, _ in series for row in values):
            axis.set_yscale('log')
        else:
            axis.set_ylim(bottom=0)
        axis.set_title(title, loc='left')
        axis.set_xlabel('Accepted step (discarded branches are separate)')
        axis.set_ylabel(unit)
        axis.grid(alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles[:4], labels[:4], loc='upper center', bbox_to_anchor=(0.5, 0.925),
                  ncol=2, frameon=False)
    figure.suptitle('T2 eight-frame training diagnostic — no independent validation', fontsize=15, y=0.985)
    figure.text(0.07, 0.02, 'Original recorded points; no smoothing, interpolation or model evaluation.', fontsize=9)
    figure.subplots_adjust(top=0.82, bottom=0.10, left=0.08, right=0.97, hspace=0.38, wspace=0.28)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, facecolor='white')
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle-dir', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path,
                        help='Plot file (.png, .pdf or .svg); default: BUNDLE/generated/t2-fit-diagnostic.png')
    parser.add_argument('--verify-only', action='store_true', help='Verify evidence without importing Matplotlib or writing files')
    args = parser.parse_args()
    require(not (args.verify_only and args.output is not None), '--output and --verify-only cannot be combined')
    root = args.bundle_dir.resolve()
    require(root.is_dir(), 'Bundle directory does not exist')
    rows, branches, report, protected_paths = verify_bundle(root)
    if not args.verify_only:
        output = (args.output or root/'generated/t2-fit-diagnostic.png').resolve()
        require(output.suffix.lower() in ('.png', '.pdf', '.svg'), 'Plot output must be PNG, PDF or SVG')
        require(output not in protected_paths, 'Plot output would overwrite verified input evidence')
        plot_curves(rows, branches, output)
        report['plot'] = str(output)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError, ImportError) as error:
        print('Verification/plotting failed: ' + str(error), file=sys.stderr)
        raise SystemExit(1)
