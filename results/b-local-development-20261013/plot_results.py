"""Plot archived aggregate scores only; no model/data inference dependencies."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
with (root/'training-curve-snapshots.csv').open(newline='',encoding='utf8') as stream:
    records=list(csv.DictReader(stream))
fig,axes=plt.subplots(1,2,figsize=(10,3.8),constrained_layout=True)
for ax,prefix,title in zip(axes,('B16','development40'),('B16 local training','Old development40 training')):
    selected={int(row['accepted_updates']):row for row in records if row['experiment'].startswith(prefix)}
    steps=sorted(selected)
    for key,label in [('force_RMSE_eV_per_A','All-atom force'),('CO2_force_RMSE_eV_per_A','CO2 force')]:
        ax.plot(steps,[float(selected[x][key]) for x in steps],marker='.',label=label)
    ax.axhline(.03,color='C0',linestyle='--',linewidth=1,label='All-atom training threshold')
    ax.axhline(.05,color='C1',linestyle=':',linewidth=1,label='CO2 group threshold (reference)')
    ax.set(xlabel='Accepted LBFGS updates',ylabel='Component RMSE (eV/Å)',title=title)
    ax.grid(alpha=.2);ax.legend(fontsize=7)
fig.suptitle('Recorded training scores only; B-H test and old40 qualification failed',fontsize=10)
fig.savefig(root/'training-snapshots.png',dpi=160)
fig.savefig(root/'training-snapshots.svg')
