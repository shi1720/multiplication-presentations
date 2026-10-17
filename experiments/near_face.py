"""Full-dimensional near-face perturbation audit; fixed seed and exact bounds.

The theoretical certificate covers every generated point; membership is checked
independently by constructing both floating symmetric matrices.
"""
from pathlib import Path
import sys,json,csv
import numpy as np
from scipy.special import gamma
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from presentations import splits,certificate_matrix
root=Path(__file__).resolve().parents[1]; rng=np.random.default_rng(6141720)
r=.25; rows=[]; N=2000
for d in [4,6,8]:
 n=d-3
 for delta in [.1,.03,.01,.003,.001]:
  certified=0; bmin=1.; cmin=1.
  for case in range(N):
   # Uniform in n-dimensional L1 ball: simplex proportions and radial CDF.
   p=rng.exponential(size=n); p/=p.sum()
   z=r*p*rng.random()**(1/n)*rng.choice([-1,1],size=n)
   w=rng.uniform(-np.sqrt(delta),np.sqrt(delta))
   u=np.concatenate(([0,.5-delta,w],z))
   ce=np.linalg.eigvalsh(certificate_matrix(u,splits(d)))[0]
   be=np.linalg.eigvalsh(certificate_matrix(u,splits(d,'balanced')))[0]
   assert be>=-1e-12
   assert (ce>=0)==(np.dot(u[2:],u[2:])<=2*delta)
   certified+=ce>=0; bmin=min(bmin,be);cmin=min(cmin,ce)
  ball=np.pi**((n+1)/2)/gamma((n+1)/2+1)*(2*delta)**((n+1)/2)
  prism=2*np.sqrt(delta)*2**n*r**n/gamma(n+1)
  rows.append(dict(d=d,delta=delta,cases=N,chain=int(certified),balanced=N,
                   chain_fraction=certified/N,chain_fraction_upper_bound=min(1,ball/prism),
                   min_balanced_eigenvalue=bmin,min_chain_eigenvalue=cmin))
with (root/'results'/'near_face.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
(root/'results'/'near_face_summary.json').write_text(json.dumps({
 'seed':6141720,'total_cases':N*len(rows),'balanced_certified':N*len(rows),
 'chain_certified':sum(x['chain'] for x in rows),'rows':rows},indent=2))
print(json.dumps(rows,indent=2))
