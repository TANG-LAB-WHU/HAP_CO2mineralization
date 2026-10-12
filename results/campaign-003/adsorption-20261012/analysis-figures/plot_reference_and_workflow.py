"""Render two development figures from saved scalar data; no model execution."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT/'data/reference-workflow.json').read_text(encoding='utf-8'))
OUT = ROOT/'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
INK='#233747'; BLUE='#0072B2'; RED='#B3482E'; GREEN='#007F6C'; PALE='#EFF4F7'

def box(ax,x,y,w,h,title,body,color=BLUE,fontsize=10):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.014,rounding_size=0.02',
        facecolor=PALE,edgecolor=color,linewidth=1.4))
    ax.text(x+.025,y+h-.06,title,ha='left',va='top',weight='bold',color=color,fontsize=fontsize+1)
    ax.text(x+.025,y+h-.14,body,ha='left',va='top',color=INK,fontsize=fontsize,linespacing=1.35)

def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=220,bbox_inches='tight')
    plt.close(fig)

def reference():
    d=DATA['A_reference'];e=d['energy_Ha'];energy=d['E_ads_A_eV_per_cell_and_per_CO2']
    fig=plt.figure(figsize=(8.4,5.6));fig.suptitle('A: matched relaxed references',x=.09,ha='left',weight='bold',color=INK,fontsize=15)
    ax=fig.add_axes([.065,.49,.90,.40]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    for x,key,title in [(.02,'A91','Adsorbed A91'),(.35,'clean88','Clean slab'),(.68,'gas3','Periodic CO2')]:
        box(ax,x,.14,.28,.79,title,f"{e[key]:.10f} Ha\nAll four optimizer\ncriteria passed",GREEN,10)
    ax.text(.327,.51,'−',ha='center',va='center',fontsize=22,color=INK)
    ax.text(.67,.51,'−',ha='center',va='center',fontsize=22,color=INK)
    ax2=fig.add_axes([.16,.33,.73,.11]);ax2.axvline(0,color=INK,lw=1)
    ax2.barh([0],[energy],height=.35,color=BLUE)
    ax2.set(xlim=(-.43,.025),yticks=[],xticks=[-.4,-.3,-.2,-.1,0],xlabel='Electronic adsorption energy (eV per CO2)')
    ax2.text(energy/2,0,f'{energy:.6f} eV',ha='center',va='center',color='white',weight='bold',fontsize=12)
    ax2.spines['left'].set_visible(False)
    fig.text(.09,.17,'Specified dry parallel-OH slab and coverage; matched cell and PBE-D3(BJ)/Z dipole.',fontsize=9.5,color=INK)
    fig.text(.09,.115,'Periodic gas limit unverified · No BSSE correction · No vibrational minimum test',fontsize=9.5,color=INK)
    fig.text(.09,.055,'B remains partial and is excluded: this figure makes no A/B ranking.',fontsize=10,color=RED,weight='bold')
    save(fig,'c003-reference-energy')

def workflow():
    d=DATA['workflow'];fig,ax=plt.subplots(figsize=(8.4,6.5));fig.subplots_adjust(left=.025,right=.985,top=.91,bottom=.075)
    ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    fig.suptitle('91-atom workflow: execution and accuracy are separate',x=.055,ha='left',fontsize=14,weight='bold',color=INK)
    box(ax,.02,.66,.29,.29,'1  CP2K reference data','Short DFT trajectories\nSame starting basin\nReference QA passed',GREEN)
    box(ax,.355,.66,.29,.29,'2  DeepMD v2: failed',f"{d['train_frames']} training frames\nForce RMSE {d['train_force_RMSE_eV_A']:.5f}\n> {d['train_force_limit_eV_A']:.2f} eV/Å",RED)
    box(ax,.69,.66,.29,.29,'3  Diagnostic: failed',f"R9/R10: {d['diagnostic_frames']} frames\nForce RMSE {d['diagnostic_force_RMSE_eV_A']:.5f}\n> {d['diagnostic_force_limit_eV_A']:.2f} eV/Å",RED)
    ax.annotate('',xy=(.347,.81),xytext=(.317,.81),arrowprops={'arrowstyle':'->','color':INK})
    ax.annotate('',xy=(.682,.81),xytext=(.652,.81),arrowprops={'arrowstyle':'->','color':RED,'linestyle':'--'})
    ax.text(.667,.625,'Separate diagnostic authorization',ha='center',fontsize=8.5,color=RED)
    box(ax,.06,.20,.40,.30,'4  LAMMPS exploration',f"{d['MD_steps']} steps × 0.25 fs = {d['MD_duration_fs']} fs\nInterface + geometry checks passed\nForce qualification still failed",GREEN)
    box(ax,.54,.20,.40,.30,'5  CP2K 10 fs probe',f"DFT reference QA passed\nCO2 force RMSE {d['probe_CO2_force_RMSE_eV_A']:.6f}\n> {d['probe_CO2_limit_eV_A']:.2f} eV/Å: screen failed",RED)
    ax.plot([.835,.835,.255],[.65,.555,.555],color=RED,linestyle='--',linewidth=1.1)
    ax.annotate('',xy=(.255,.51),xytext=(.255,.555),arrowprops={'arrowstyle':'->','color':RED,'linestyle':'--'})
    ax.text(.49,.58,'Separately authorized exploration; no model qualification',ha='center',fontsize=8.3,color=RED)
    ax.annotate('',xy=(.53,.35),xytext=(.47,.35),arrowprops={'arrowstyle':'->','color':INK})
    ax.text(.74,.14,'20 fs DFT probe skipped by the frozen rule',ha='center',fontsize=9,color=RED,weight='bold')
    ax.text(.5,.05,'Observed campaign-002 v2 route. R4 remains sealed; no production potential.',ha='center',fontsize=9.3,color=INK)
    save(fig,'c003-validation-workflow')

if __name__=='__main__':
    reference();workflow()
    print('Rendered c003-reference-energy and c003-validation-workflow from compact saved data.')
