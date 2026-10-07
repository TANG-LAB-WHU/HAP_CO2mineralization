"""Render the actual successful short window, with its longer-time failure explicit."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from matplotlib.lines import Line2D
from ase.io import read
from plot_refinement import COLORS
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run;d=r/'lammps-recovered'
    if (r/'reproducibility-manifest.json').exists():raise RuntimeError('Captured recovery is immutable')
    audit=json.loads((r/'recovered-dynamics.json').read_text());assert audit['geometry_gate'] and audit['saved_frames']==101
    path=d/'trajectory.extxyz';frames=read(path,index=':');t=np.load(d/'thermo.npy');ids=np.arange(0,101,2)
    assert len(frames)==len(t)==101 and abs(t[-1,1]-.25)<1e-12
    positions=np.concatenate([frames[i].positions for i in ids]);low=positions.min(axis=0)-1;high=positions.max(axis=0)+1
    fig,axes=plt.subplots(1,2,figsize=(10,6));fig.subplots_adjust(top=.79,bottom=.15);artists=[]
    colors=[COLORS[s] for s in frames[0].get_chemical_symbols()]
    for axis,view in zip(axes,[(0,2),(0,1)]):
        artists.append(axis.scatter([],[],s=50,edgecolors='white',linewidths=.3))
        axis.set(xlim=(low[view[0]],high[view[0]]),ylim=(low[view[1]],high[view[1]]),xlabel='xyz'[view[0]]+' (A)',ylabel='xyz'[view[1]]+' (A)');axis.set_aspect('equal')
    title=fig.suptitle('',fontsize=13)
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=s) for s,c in COLORS.items()],ncol=5,loc='upper center',bbox_to_anchor=(.5,.90))
    fig.text(.025,.035,'Actual 250 fs demonstration only. This same model fails near 395 fs in the separate 1 ps run.',fontsize=9,color='#9c6416')
    def update(k):
        i=ids[k]
        for artist,view in zip(artists,[(0,2),(0,1)]):artist.set_offsets(frames[i].positions[:,view]);artist.set_color(colors)
        title.set_text(f'STABLE SHORT WINDOW ONLY | t = {t[i,1]*1000:.0f} fs | T = {t[i,2]:.1f} K\nNVE; dt = 0.25 fs; all atoms mobile; exploratory force gate still fails')
        return artists+[title]
    movie=FuncAnimation(fig,update,frames=len(ids),interval=125);movie.save(r/'13_stable250fs_demo.gif',writer=PillowWriter(fps=8));plt.close(fig)
    record={'frames':len(ids),'source_frame_indices':ids.tolist(),'source_time_fs':(t[ids,1]*1000).tolist(),'interpolation':False,
        'scope':'Actual separate 250 fs run; no long-time stability claim; same model fails 1 ps geometry and retrospective force gates.',
        'sources':{str(p.relative_to(r)):sha(p) for p in [path,d/'thermo.npy',r/'recovered-dynamics.json']},'script_sha256':sha(Path(__file__)),'manuscript_eligible':False}
    (r/'stable-demo-manifest.json').write_text(json.dumps(record,indent=2)+'\n');print('Rendered 51 actual stable-window frames with explicit later-failure label')
if __name__=='__main__':main()
