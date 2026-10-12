"""Verify saved hashes and endpoint arithmetic; never evaluates a model."""
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
def read_json(name): return json.loads((ROOT / name).read_text(encoding='utf8'))
def read_csv(name):
    with (ROOT / name).open(encoding='utf8', newline='') as stream:
        return list(csv.DictReader(stream))
def close(actual, expected, name):
    if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError(f'{name}: {actual} != {expected}')

manifest = read_json('manifest.json')
for name, entry in manifest['files'].items():
    path = (ROOT / name).resolve()
    if not path.is_relative_to(ROOT) or path.is_symlink():
        raise ValueError('Unsafe package path: ' + name)
    raw = path.read_bytes()
    if len(raw) != entry['bytes'] or hashlib.sha256(raw).hexdigest() != entry['sha256']:
        raise ValueError('File changed: ' + name)

summary = read_json('summary.json')
rows = read_csv('data/endpoint-10fs-forces.csv')
assert len(rows) == 91 and [int(r['atom_id']) for r in rows] == list(range(1, 92))
assert [r['species'] for r in rows[-3:]] == ['C', 'O', 'O']
errors = []
for row in rows:
    vector = []
    for axis in ('x', 'y', 'z'):
        value = float(row[f'ML_{axis}_eV_A']) - float(row[f'DFT_{axis}_eV_A'])
        close(value, float(row[f'error_ML_minus_DFT_{axis}_eV_A']), 'raw force difference')
        vector.append(value)
    errors.append(vector)
flat = [x for row in errors for x in row]
metrics = summary['endpoint_10fs']['metrics']
close(math.sqrt(sum(x*x for x in flat) / 273), metrics['force_RMSE_eV_A'], 'force RMSE')
close(max(abs(x) for x in flat), metrics['max_abs_force_component_error_eV_A'], 'maximum force error')
co2 = [x for row in errors[-3:] for x in row]
close(math.sqrt(sum(x*x for x in co2) / 9), metrics['CO2_force_RMSE_eV_A'], 'CO2 force RMSE')
energy = read_json('data/endpoint-10fs-energy.json')
close(abs(energy['ML_potential_energy_eV'] - energy['energy_eV']) / 91, metrics['energy_abs_error_eV_atom'], 'energy error')
thermo = read_csv('data/thermo.csv')
assert len(thermo) == 81 and [int(float(r['step'])) for r in thermo] == list(range(81))
initial = float(thermo[0]['total_eV'])
close(max(abs(float(r['total_eV']) - initial) for r in thermo), summary['MD']['energy_max_abs_change_eV'], 'energy departure')
assert summary['endpoint_20fs']['DFT_result'] is None
assert summary['endpoint_10fs']['original_endpoint_gates_passed'] is False
assert summary['endpoint_10fs']['exploratory_gates']['CO2_force_RMSE'] is False
assert summary['manuscript_eligible'] is False and summary['production_released'] is False
assert summary['R4_numeric_read_or_hashed'] is False
print(f'Verified {len(manifest["files"])} package files, 273 uncorrected endpoint force components and 81 thermo frames. Original force failure retained. No model evaluation or R4 access.')
