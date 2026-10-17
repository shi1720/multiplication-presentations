from pathlib import Path
import sys,json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from presentations import *
ROOT=Path(__file__).resolve().parents[1]; F=ROOT/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
 'axes.spines.right':False,'pdf.fonttype':42,'figure.dpi':150,'savefig.bbox':'tight'})
blue='#145DA0'; orange='#C65E25'; green='#24866B'; gray='#6C7480'
def save(name):
 plt.savefig(F/(name+'.pdf')); plt.savefig(F/(name+'.png'),dpi=180); plt.close()

fig,axes=plt.subplots(1,2,figsize=(10.2,3.65))
w=np.linspace(-3,3,401); b=np.array([w_boundary(4,x)[0] for x in w])
axes[0].fill_between(w,(1-w*w)/2,b,color=orange,alpha=.16,label='valid but chain not certified')
axes[0].plot(w,b,color=blue,lw=2.2,label='true cell / full closure')
axes[0].plot(w,(1-w*w)/2,color=orange,ls='--',label='chain and balanced (z = 0)')
axes[0].set(xlabel=r'$w=u_3$ (with $u_4=0$)',ylabel=r'largest certifiable $u_2$',title='A  |  Curved normal slice')
axes[0].legend(frameon=False,fontsize=8,loc='lower center')
z=np.linspace(-3,1.5,401); b=np.array([normal_boundary([0,x])[0] for x in z])
axes[1].plot(z,b,color=blue,lw=2.2,label='true Voronoi boundary')
axes[1].plot(z,(1-z*z)/2,color=orange,ls='--',label='chain')
axes[1].plot(z[z<=.5],np.full(sum(z<=.5),.5),color=green,ls=':',lw=3,label='balanced')
axes[1].scatter([-2,1],[.5,-1],color=[green,orange],zorder=5,s=28)
axes[1].axvline(.5,color=gray,alpha=.5,lw=.8)
axes[1].text(.58,-3.2,'balanced has no\ncertificate for z > 1/2',fontsize=8)
axes[1].set(xlabel=r'$z=u_4$ (with $u_3=0$)',title='B  |  An unbounded exposed face')
axes[1].legend(frameon=False,fontsize=8,loc='lower left')
fig.tight_layout(w_pad=2);save('normal_slices')

fig,axes=plt.subplots(1,2,figsize=(10.2,3.4))
for d in [3,4,6,8]:
 w=np.geomspace(.008,.2,80); gap=np.array([w_boundary(d,x)[0]-(1-x*x)/2 for x in w])
 axes[0].loglog(w,gap,label=f'd = {d}')
w=np.geomspace(.008,.2,80)
axes[0].loglog(w,.5*w**4,color='black',ls=':',label=r'$w^4/2$')
axes[0].set(xlabel=r'$|u_3|$',ylabel='vertical certificate gap',title='A  |  Fourth-order contact')
axes[0].legend(frameon=False,fontsize=8)
d=np.arange(3,13)
axes[1].plot(d,d-3,'o-',color=blue,label='true face and balanced presentation')
axes[1].plot(d,np.zeros(len(d)),'s--',color=orange,label='chain presentation')
axes[1].set(xlabel='ambient dimension d',ylabel='exposed-face dimension',title='B  |  Dimension lost at fixed matrix size',xticks=d[::2])
axes[1].legend(frameon=False,fontsize=8)
fig.tight_layout(w_pad=2);save('contact_and_dimension')

fig,axes=plt.subplots(1,3,figsize=(10.2,3.5),sharey=True)
from matplotlib.colors import ListedColormap
colors=['#ECEFF2','#B7D7EB','#E0A47A','#68AD96','#384C61']
labels=['outside cell','valid; neither','chain only','balanced only','both']
rows=list(csv.DictReader(open(ROOT/'results'/'rational_grid.csv')))
from fractions import Fraction
for ax,w0 in zip(axes,[0,.5,1]):
 sub=[r for r in rows if float(Fraction(r['w']))==w0]
 for r in sub:
  T=r['true_cell']=='True'; C=r['chain']=='True'; B=r['balanced']=='True'
  cls=0 if not T else (4 if C and B else (2 if C else (3 if B else 1)))
  ax.scatter(float(Fraction(r['z'])),float(Fraction(r['v'])),s=55,c=[colors[cls]],marker='s')
 ax.set(xlabel=r'$u_4$',title=rf'$u_3={w0:g}$',xticks=[-2,-1,0,1,2])
axes[0].set_ylabel(r'$u_2$')
handles=[plt.Line2D([],[],marker='s',ls='',color=c,label=l) for c,l in zip(colors,labels)]
fig.legend(handles=handles,loc='lower center',ncol=5,frameon=False,bbox_to_anchor=(.5,-.05),fontsize=8)
fig.tight_layout();save('exact_grid')

fig,ax=plt.subplots(figsize=(7.4,3))
z=np.linspace(-4,0,41);delta=np.array([.001,.01,.05,.1,.25])
A=np.array([[int(dd>=zz*zz/2-1e-12) for zz in z] for dd in delta])
ax.imshow(A,origin='lower',aspect='auto',extent=[-4.05,.05,-.5,4.5],cmap=ListedColormap(['#E8EEF4',blue]),vmin=0,vmax=1)
ax.set_yticks(range(5),[str(x) for x in delta]);ax.set(xlabel=r'fourth-channel perturbation $u_4$',ylabel=r'offset $delta=1/2-u_2$')
ax.set_title('Chain certifies 20/205; balanced certifies all 205 synthetic inputs',fontsize=11)
fig.tight_layout();save('synthetic_calibration')
print('Wrote four vector figures and PNG previews')
