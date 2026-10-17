"""Deterministic exact audit, numerical SDP comparison, and synthetic use case."""
import sys,json,csv,time,platform,warnings
from pathlib import Path
import numpy as np
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from presentations import *
ROOT=Path(__file__).resolve().parents[1]; R=ROOT/'results'
def save(name,rows):
 with (R/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

# Exact rational grid: no optimizer is used to decide membership.
rows=[]
for v in [-2,-1,sp.Rational(-1,2),0,sp.Rational(1,4),sp.Rational(49,100),sp.Rational(1,2)]:
 for w in [-1,sp.Rational(-1,2),0,sp.Rational(1,2),1]:
  for z in [-2,-1,sp.Rational(-1,2),0,sp.Rational(1,4),sp.Rational(1,2),sp.Rational(3,4),1,2]:
   u=[0,v,w,z]; true=exact_nonnegative(u)
   chain=rational_psd(exact_matrix(u,splits(4)))
   balanced=rational_psd(exact_matrix(u,splits(4,'balanced')))
   assert not(chain or balanced) or true
   rows.append(dict(v=str(v),w=str(w),z=str(z),true_cell=true,chain=chain,balanced=balanced))
save('rational_grid.csv',rows)
grid={'total':len(rows),'true_cell':sum(x['true_cell'] for x in rows),
      'chain':sum(x['chain'] for x in rows),'balanced':sum(x['balanced'] for x in rows),
      'union':sum(x['chain'] or x['balanced'] for x in rows),
      'chain_only':sum(x['chain'] and not x['balanced'] for x in rows),
      'balanced_only':sum(x['balanced'] and not x['chain'] for x in rows),
      'false_positives':0}
print('exact grid',grid,flush=True)

# Numerical closure checks: 12 fixed slice points and 12 random tails per d.
rng=np.random.default_rng(1720406); sdprows=[]
for d in range(3,9):
 for case in range(24):
  tail=np.zeros(d-2)
  if case<12: tail[0]=np.linspace(-3,3,12)[case]
  else: tail=rng.uniform(-.6,.6,d-2)
  val,t=normal_boundary(tail)
  sol=closure_boundary(tail)
  row={'d':d,'case':case,'tail':json.dumps(tail.tolist()),'scalar_boundary':val,
       'parameter':t,**sol,'absolute_error':abs(sol['boundary']-val)}
  assert row['absolute_error']<2e-6,row
  assert row['coefficient_residual']<2e-7,row
  assert row['min_eigenvalue']>-2e-7,row
  sdprows.append(row)
 print('closure d',d,'done',flush=True)
save('closure_checks.csv',sdprows)

# Synthetic normal perturbations of a four-channel monomial sensor.
app=[]
for delta in [.001,.01,.05,.1,.25]:
 for z in np.linspace(-4,0,41):
  u=[0,.5-delta,0,z]
  app.append({'delta':delta,'z':z,'chain':bool(np.linalg.eigvalsh(certificate_matrix(u,splits(4)))[0]>=-1e-12),
              'balanced':bool(np.linalg.eigvalsh(certificate_matrix(u,splits(4,'balanced')))[0]>=-1e-12)})
save('synthetic_calibration.csv',app)

# Analytic dimension-dependent boundary samples, stable root calculation.
bounds=[]
for d in range(3,9):
 for w in np.concatenate(([0],np.geomspace(.002,10,160))):
  b,s=w_boundary(d,w)
  bounds.append({'d':d,'w':w,'true_boundary':b,'chain_boundary':(1-w*w)/2,'parameter':s})
save('boundaries.csv',bounds)

import cvxpy,scipy
summary={'seed':1720406,'python':platform.python_version(),'platform':platform.platform(),
         'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,'cvxpy':cvxpy.__version__,
         'exact_grid':grid,'sdp_comparisons':len(sdprows),
         'max_sdp_error':max(x['absolute_error'] for x in sdprows),
         'max_coefficient_residual':max(x['coefficient_residual'] for x in sdprows),
         'min_sdp_eigenvalue':min(x['min_eigenvalue'] for x in sdprows),
         'optimal_inaccurate_count':sum(x['status']!='optimal' for x in sdprows),
         'synthetic_cases':len(app),'synthetic_chain':sum(x['chain'] for x in app),
         'synthetic_balanced':sum(x['balanced'] for x in app)}
(R/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
