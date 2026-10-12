"""Replay saved report arithmetic; no scientific execution or raw/sealed data access."""
from pathlib import Path
import csv
import hashlib
import json
import math

root = Path(__file__).resolve().parent
load = lambda name: json.loads((root / name).read_text(encoding='utf8'))
manifest = load('manifest.json')
for name, record in manifest['files'].items():
    path = (root / name).resolve()
    if not path.is_relative_to(root) or path.is_symlink():
        raise ValueError('Invalid package path: ' + name)
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != record['sha256'] or len(raw) != record['bytes']:
        raise ValueError('Changed package file: ' + name)

ads = load('data/A-conditional-adsorption-energy.json')
energy = ads['energy_Ha']
value = (energy['A91'] - energy['clean88'] - energy['gas3']) * 27.211386245988
assert math.isclose(value, ads['E_ads_A_eV_per_cell_and_per_CO2'], rel_tol=0, abs_tol=1e-12)
assert ads['B_used'] is False and ads['manuscript_eligible'] is False
endpoints = load('endpoints-summary.json')
for name in ('A91', 'clean88', 'gas3'):
    assert endpoints[name]['optimization_converged'] is True
    assert endpoints[name]['final_energy_only'][-1]['energy_Ha'] == energy[name]
    assert endpoints[name]['final_energy_only'][-1]['excluded_from_force_labels'] is True
b = endpoints['B91-continuation1']
assert b['status'] == 'partial_budget' and b['optimization_converged'] is False
assert b['endpoint_eligible_for_ranking'] is False
assert b['geometry_updates_with_optimizer_output'] == 20
with (root / 'data/optimization-snapshots.csv').open(encoding='utf8', newline='') as stream:
    rows = list(csv.DictReader(stream))
for name, item in endpoints.items():
    selected = [r for r in rows if r['segment'] == name]
    assert len(selected) == item['force_evaluation_count']
    assert sum(r['complete'] == 'True' for r in selected) == item['complete_force_evaluation_count']
    for row in selected:
        if row['complete'] == 'False':
            assert row['energy_Ha'] == '' and row['maximum_force_Ha_per_bohr'] == ''
assert load('omitted-artifacts.json')['R4_numeric_files_opened_or_hashed'] is False
print(f'Verified {len(manifest["files"])} files, Eads(A)={value:.14f} eV, converged A/clean/gas, and partial unranked B. No new computation.')
