"""Render measured development tables only; no model or scientific execution."""
from pathlib import Path
import csv, json, os, tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'hap91-closeout-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root=Path(__file__).resolve().parent
out=root/'figures';out.mkdir(exist_ok=True)
def read(n):
    with (root/'data'/n).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
scores=read('force-and-energy.csv');elements=read('endpoint-elements.csv')
s=json.loads((root/'summary.json').read_text(encoding='utf-8'))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(1,3,figsize=(13,4.8))
groups=[scores[0],scores[1],scores[-1]]
x=np.arange(3); vals=[float(t['force_RMSE_eV_A']) for t in groups]
ax[0].bar(x,vals,color=['#b04a40']*3,width=.57)
ax[0].scatter(x,[.03,.05,.05],marker='_',s=900,c='black',label='Original gate',zorder=5)
ax[0].set_xticks(x,['Train\n90 frames','R9/R10\n16 frames','10-fs probe\n1 endpoint'])
ax[0].set_ylim(0,.12);ax[0].set_ylabel('Force-component RMSE (eV/Å)')
ax[0].set_title('A  Original force gates remain unmet',loc='left',fontsize=11)
for i,v in enumerate(vals):ax[0].text(i,v+.003,f'{v:.5f}',ha='center',fontsize=9)
ax[0].legend(frameon=False,fontsize=9,loc='upper left')
co=[float(t['CO2_force_RMSE_eV_A']) for t in groups[1:]]
ax[1].bar([0,1],co,color=['#61788b','#b04a40'],width=.57)
ax[1].axhline(.05,color='black',ls='--',lw=1,label='Endpoint exploratory limit')
ax[1].set_xticks([0,1],['R9/R10\ndiagnostic','10-fs\nDFT probe'])
ax[1].set_ylim(0,.08);ax[1].set_ylabel('CO₂ force-component RMSE (eV/Å)')
ax[1].set_title('B  CO₂ endpoint screen stops the run',loc='left',fontsize=11)
for i,v in enumerate(co):ax[1].text(i,v+.002,f'{v:.5f}',ha='center',fontsize=9)
ax[1].legend(frameon=False,fontsize=8,loc='upper left')
xx=np.arange(5);ev=[float(t['force_RMSE_eV_A']) for t in elements]
ax[2].bar(xx,ev,color=['#b04a40' if v>.075 else '#61788b' for v in ev],width=.65)
ax[2].axhline(.075,color='black',ls='--',lw=1,label='Original 0.075')
ax[2].axhline(.15,color='#61788b',ls=':',lw=1,label='Exploratory 0.15')
ax[2].set_xticks(xx,[t['element'] for t in elements]);ax[2].set_ylim(0,.18)
ax[2].set_ylabel('10-fs element force RMSE (eV/Å)')
ax[2].set_title('C  Substrate errors also remain',loc='left',fontsize=11)
ax[2].legend(frameon=False,fontsize=8,loc='upper right')
fig.suptitle('HAP91 coverage-v2: completed workflow, unqualified potential',fontsize=15,y=1.02)
fig.text(.5,-.03,'DEVELOPMENT DIAGNOSTIC • No uncertainty estimate • Different sample sets; not a controlled before/after comparison\n20-fs DFT skipped by the frozen rule; no value is imputed. All original failures retained.',ha='center',fontsize=9,color='#555555')
fig.tight_layout()
for ext in ('png',):fig.savefig(out/f'force-gates.{ext}',dpi=180,bbox_inches='tight',metadata={'Title':'HAP91 development force diagnostic'})
plt.close(fig)
t=read('thermo.csv');times=np.array([float(r['time_ps'])*1000 for r in t]);E=np.array([float(r['total_eV']) for r in t]);T=np.array([float(r['temperature_K']) for r in t])
assert len(t)==81 and float(t[-1]['step'])==80
fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].plot(times,(E-E[0])*1000/91,color='#286c8e');ax[0].set_ylabel('Total energy − initial (meV/atom)')
ax[0].set_title('A  Small total-energy departure over 20 fs',loc='left',fontsize=11)
ax[0].set_ylim(-.0003,.0085)
ax[0].text(.03,.94,'Maximum |ΔE| = 0.005717 meV/atom\nOriginal 1 meV/atom limit is outside this zoom',transform=ax[0].transAxes,va='top',fontsize=9)
ax[1].plot(times,T,color='#a05d27');ax[1].set_ylabel('Instantaneous temperature (K)')
ax[1].set_title('B  300 K initialization is not equilibrium',loc='left',fontsize=11)
for a in ax:
    a.set_xlabel('Time (fs)');a.axvline(10,color='#b04a40',ls=':',lw=1);a.set_xlim(0,20)
fig.suptitle('HAP91 NVE diagnostic • one seed, 80 steps, 81 correlated frames',fontsize=13,y=1.02)
fig.text(.5,-.03,'DEVELOPMENT DIAGNOSTIC • Stable integration does not establish DFT force accuracy\nRed line: the only DFT-checked endpoint; its CO₂ force screen failed.',ha='center',fontsize=9,color='#555555')
fig.tight_layout()
for ext in ('png',):fig.savefig(out/f'nve-diagnostic.{ext}',dpi=180,bbox_inches='tight',metadata={'Title':'HAP91 development integration diagnostic'})
plt.close(fig)
print('Rendered two measured diagnostic figures as PNG; no model evaluation.')
