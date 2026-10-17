"""Exact weighted identities and PSD constructions, 400 deterministic cases."""
import sys,json
from pathlib import Path
import numpy as np
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from presentations import all_presentations,exact_matrix,rational_psd
root=Path(__file__).resolve().parents[1];rng=np.random.default_rng(614)
t=sp.Symbol('t');cases=0
for d in range(3,11):
 options=list(all_presentations(d))
 for rep in range(50):
  P=options[int(rng.integers(len(options)))];b=[sp.Rational(int(x),3) for x in rng.integers(1,20,size=d)]
  u=[sp.Integer(0),b[0]/(2*b[1])]+[
    sp.Rational((-1)**k,100*d)*min(b)/b[k-1] if P[k-2][0]>1 else sp.Integer(0)
    for k in range(3,d+1)]
  H=sp.diag(*b)
  for k,(i,j) in enumerate(P,start=2):
   if i==j:H[i-1,i-1]-=2*b[k-1]*u[k-1]
   else:
    H[i-1,j-1]-=b[k-1]*u[k-1];H[j-1,i-1]-=b[k-1]*u[k-1]
  z=sp.Matrix([t**k for k in range(1,d+1)])
  direct=sum(b[k-1]*(t**(2*k)-2*u[k-1]*t**k) for k in range(1,d+1))
  assert sp.expand((z.T*H*z)[0]-direct)==0
  assert rational_psd(H)
  for k,(i,j) in enumerate(P,start=2):
   if k>=3 and i==1:
    M=H.copy();M[0,j-1]=-sp.Rational(1,10);M[j-1,0]=M[0,j-1]
    assert not rational_psd(M)
  cases+=1
(root/'results'/'weighted_verification.json').write_text(json.dumps({'cases':cases,'seed':614,'passed':True},indent=2))
print(cases,'weighted cases passed')
