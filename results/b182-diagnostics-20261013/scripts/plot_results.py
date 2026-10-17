#!/usr/bin/env python3
"""Verify and redraw the archived B182 tables. No models or scientific engines.

Python 3.9+; Matplotlib is needed only for plotting. Paths are relative to this
script, so the bundle can be moved. No smoothing or interpolation is performed.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def col(rows, name):
    return [float(row[name]) for row in rows]


def verify(root):
    manifest = json.loads((root / 'bundle-manifest.json').read_text(encoding='utf-8'))
    for name, record in manifest['files'].items():
        data = (root / name).read_bytes()
        if hashlib.sha256(data).hexdigest() != record['sha256'] or len(data) != record['bytes']:
            raise ValueError('Bundled file differs: ' + name)
    for name, count in [('dft4fs-curve.csv', 17), ('ml85fs-curve.csv', 341), ('ml85fs-OH-distances.csv', 341)]:
        rows = read_csv(root / 'tables' / name)
        if len(rows) != count or col(rows, 'step') != list(range(count)):
            raise ValueError('Incomplete or reordered table: ' + name)
        if col(rows, 'time_fs') != [i * 0.25 for i in range(count)]:
            raise ValueError('Time/step mismatch: ' + name)
    rows = read_csv(root / 'tables/step340-all-atom-forces.csv')
    if len(rows) != 182 or col(rows, 'atom_id') != list(range(1, 183)):
        raise ValueError('Force table does not retain all 182 ordered atoms')
    return {'verified_files': len(manifest['files']), 'DFT4fs_frames': 17,
            'ML85fs_frames': 341, 'step340_force_rows': 182,
            'model_loaded': False, 'scientific_calculation_performed': False}


def plot(root, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .22,
                         'savefig.dpi': 160, 'figure.facecolor': 'white'})
    output.mkdir(parents=True, exist_ok=True)
    produced = []

    def save(fig, filename, subtitle):
        fig.text(.5, .013, subtitle, ha='center', va='bottom', fontsize=9, color='#7b2020')
        fig.savefig(output / filename, bbox_inches='tight', metadata={'Software': 'B182 archived-table plotter'})
        plt.close(fig)
        produced.append(filename)

    d = read_csv(root / 'tables/dft4fs-curve.csv')
    t = col(d, 'time_fs')
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), layout='constrained')
    fig.suptitle('Independent direct CP2K trajectory: 182 atoms, 4 fs, NVE')
    ax[0,0].plot(t, col(d, 'temperature_K'), marker='o', ms=3, color='#2166ac')
    ax[0,0].set_ylabel('Temperature (K)')
    ax[0,1].plot(t, col(d, 'total_energy_change_meV_per_atom'), marker='o', ms=3, color='#7b3294')
    ax[0,1].set_ylabel('Signed total-energy change (meV/atom)')
    for key in [k for k in d[0] if k.startswith('H') and k.endswith('_A')]:
        ax[1,0].plot(t, col(d, key), marker='.', ms=3, label=key[:-2].replace('_', '-'))
    ax[1,0].set_ylabel('Original O-H distance (Angstrom)')
    ax[1,0].legend(fontsize=7, ncol=2)
    for i, color in [(1,'#2166ac'),(2,'#d6604d')]:
        for j, style in [(1,'-'),(2,'--')]:
            ax[1,1].plot(t,col(d,f'CO2_{i}_C_O{j}_A'),style,color=color,label=f'CO2 {i}, C-O{j}')
    ax[1,1].set_ylabel('C-O distance (Angstrom)')
    ax[1,1].legend(fontsize=8)
    for a in ax.flat:a.set_xlabel('Time (fs)')
    fig.get_layout_engine().set(rect=(0,.05,1,.96))
    save(fig,'dft4fs.png','17 recorded DFT frames; short response only. This is not the ML trajectory or equilibrium sampling.')

    m = read_csv(root / 'tables/ml85fs-curve.csv')
    t = col(m, 'time_fs')
    fig, ax = plt.subplots(3, 2, figsize=(12, 10), layout='constrained')
    fig.suptitle('Unvalidated exploratory ML, model failed validation\n182 atoms: stopped at 85 fs after O-H elongation')
    ax[0,0].plot(t,col(m,'temperature_K'),color='#2166ac')
    ax[0,0].set_ylabel('Temperature (K)')
    ax[0,1].plot(t,col(m,'signed_energy_departure_meV_per_atom'),color='#7b3294')
    ax[0,1].set_ylabel('Signed total-energy change (meV/atom)')
    for i,color in [(1,'#2166ac'),(2,'#d6604d')]:
        for j,style in [(1,'-'),(2,'--')]:
            ax[1,0].plot(t,col(m,f'CO2_{i}_C_O{j}_A'),style,color=color,label=f'CO2 {i}, C-O{j}')
        ax[1,1].plot(t,col(m,f'CO2_{i}_O_C_O_degrees'),color=color,label=f'CO2 {i}')
        ax[2,0].plot(t,col(m,f'CO2_{i}_COM_height_local_Ca_A'),color=color,label=f'CO2 {i}')
        ax[2,1].plot(t,col(m,f'CO2_{i}_Ca_O_contact_pairs_lt_3_2_A'),color=color,label=f'CO2 {i}')
    ax[1,0].set_ylabel('C-O distance (Angstrom)')
    ax[1,1].set_ylabel('O-C-O angle (degrees)')
    ax[2,0].set_ylabel('CO2 COM height above fixed local Ca (Angstrom)')
    ax[2,1].set_ylabel('Ca-O contact pairs below 3.2 Angstrom')
    for a in ax.flat:
        a.set_xlabel('Time (fs)')
        a.axvline(82.75,color='#b2182b',ls=':',lw=1)
    for a in [ax[1,0],ax[1,1],ax[2,0],ax[2,1]]:a.legend(fontsize=8)
    fig.get_layout_engine().set(rect=(0,.05,1,.95))
    save(fig,'ml85fs.png','341 recorded frames; first O-H alarm at 82.75 fs. Low energy drift does not establish DFT force accuracy.')

    o = read_csv(root / 'tables/ml85fs-OH-distances.csv')
    fig, ax = plt.subplots(figsize=(11, 6),layout='constrained')
    fig.suptitle('Unvalidated exploratory ML: original O-H distances')
    for key in [k for k in o[0] if k.startswith('H')]:
        hot = key == 'H109_O135_A'
        ax.plot(col(o,'time_fs'),col(o,key),label=key[:-2].replace('_','-'),
                lw=2.4 if hot else 1,alpha=1 if hot else .6,color='#b2182b' if hot else None)
    ax.axhline(1.4,color='black',ls='--',label='1.4 Angstrom stop threshold')
    ax.axvline(82.75,color='#b2182b',ls=':',label='First alarm: 82.75 fs')
    ax.scatter([85],[1.4986458420655033],color='#b2182b',zorder=4)
    ax.annotate('85 fs stop\nH109-O135 = 1.498646 Angstrom',xy=(85,1.4986458420655033),
                xytext=(35,1.31),arrowprops={'arrowstyle':'->'},fontsize=10)
    ax.set(xlabel='Time (fs)',ylabel='Original O-H distance (Angstrom)',xlim=(0,87))
    ax.legend(loc='upper left',fontsize=8,ncol=2)
    fig.get_layout_engine().set(rect=(0,.065,1,.95))
    save(fig,'OH-diagnostic.png','Stopped at 85 fs; model failed validation. Elongation is a trajectory diagnostic, not proof of reaction.')

    rows=read_csv(root/'tables/step340-force-metrics.csv')
    labels={'all_atoms':'All 182','CO2_1':'CO2 1','CO2_2':'CO2 2','CO2_pooled':'Both CO2',
            'OH_H109_O135':'H109 + O135','local_union_3_5A_including_OH':'Local 20','local_environment_excluding_OH':'Local 18 (no pair)'}
    f=json.loads((root/'step340-force-analysis.json').read_text())
    fig,ax=plt.subplots(1,2,figsize=(13,6))
    fig.suptitle('85 fs endpoint: full 182-atom DFT force diagnostic')
    y=list(range(len(rows)))
    ax[0].barh([v-.17 for v in y],col(rows,'force_component_RMSE_eV_A'),height=.34,label='Component RMSE',color='#2166ac')
    ax[0].barh([v+.17 for v in y],col(rows,'max_abs_component_eV_A'),height=.34,label='Maximum |component error|',color='#d6604d')
    ax[0].set_yticks(y,[labels[r['group']]for r in rows]);ax[0].invert_yaxis()
    ax[0].set_xlabel('ML minus DFT error (eV / Angstrom)');ax[0].legend(fontsize=8)
    p=f['OH_axis']['projections']
    ax[1].bar(['ML','DFT'],[p[k]['relative_acceleration_on_axis_A_fs2']for k in ['ML','DFT']],color=['#2166ac','#d6604d'])
    ax[1].axhline(0,color='black',lw=.8)
    ax[1].set_ylabel('Force-only relative acceleration on O135 -> H109\n(Angstrom / fs^2)')
    ax[1].text(.5,.98,'Both projections are restoring (negative)',ha='center',va='top',
               transform=ax[1].transAxes,fontsize=10,bbox={'facecolor':'white','alpha':.9,'edgecolor':'none'})
    fig.subplots_adjust(left=.14,right=.98,bottom=.17,top=.87,wspace=.32)
    save(fig,'step340-force-errors.png','Unvalidated ML. These are total local forces; no isolated-bond or reaction claim. One snapshot is not validation.')
    return produced


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--verify-only',action='store_true')
    args=parser.parse_args()
    result=verify(args.data_root)
    if not args.verify_only:
        result['figures']=plot(args.data_root,args.output_dir or args.data_root/'figures')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
