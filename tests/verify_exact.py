"""Exhaustive symbolic checks over all multiplication presentations, d=3..10."""
import sys,json,time
from pathlib import Path
from collections import Counter
import numpy as np
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from presentations import *
root=Path(__file__).resolve().parents[1]
start=time.perf_counter(); rows=[]; count=0; t=sp.Symbol('t')
for d in range(3,11):
    dims=Counter(); checked=0
    # Two rational points per presentation: one on its relative face interior,
    # one with every tail coordinate nonzero to detect each forbidden edge.
    for P in all_presentations(d):
        dim=face_dimension(P); dims[dim]+=1
        u=[sp.Integer(0),sp.Rational(1,2)]+[
            sp.Rational((-1)**k,8*d) if P[k-2][0]>1 else sp.Integer(0)
            for k in range(3,d+1)]
        H=exact_matrix(u,P); z=sp.Matrix([t**k for k in range(1,d+1)])
        direct=sum(t**(2*k)-2*u[k-1]*t**k for k in range(1,d+1))
        assert sp.expand((z.T*H*z)[0]-direct)==0
        assert rational_psd(H)
        assert np.linalg.eigvalsh(np.array(H).astype(float))[0]>=-1e-12
        # Every leaf-one coefficient is structurally forced to zero by H11=0.
        for k,(i,j) in enumerate(P,start=2):
            if k>=3 and i==1:
                assert H[0,j-1]==0
                test=u.copy(); test[k-1]=sp.Rational(1,8*d)
                M=exact_matrix(test,P)
                assert M.extract([0,j-1],[0,j-1]).det()<0
                assert not rational_psd(M)
        checked+=1
    assert set(dims)==set(range(d-2))
    zvar=sp.Symbol('z')
    generating=sp.Poly(sp.prod(1+(k//2-1)*zvar for k in range(3,d+1)),zvar)
    assert all(generating.nth(j)==dims[j] for j in range(d-2))
    rows.append({'d':d,'presentations':checked,'face_dimension_counts':dict(sorted(dims.items()))})
    count+=checked
# Independently establish exact nonnegativity and incomparable sparse examples.
examples=[([0,sp.Rational(1,2),0,-2],False,True),
          ([0,-1,0,1],True,False),
          ([0,sp.Rational(1,2),0,sp.Rational(3,4)],False,False)]
for u,c,b in examples:
    assert rational_psd(exact_matrix(u,splits(4)))==c
    assert rational_psd(exact_matrix(u,splits(4,'balanced')))==b
assert exact_nonnegative(examples[0][0])
assert exact_nonnegative(examples[1][0])
assert not exact_nonnegative(examples[2][0])
# Constant, double-real-root, odd-real-root, and no-real-root checks of oracle.
for u,answer in [([0,0,0],True),([0,sp.Rational(1,2),0],True),
                 ([0,1,0],False),([0,0,1],True),
                 ([0,sp.Rational(1,2),1],False)]:
    assert exact_nonnegative(u)==answer
out={'total_presentations':count,'symbolic_identity_checks':count,
     'rows':rows,'seconds':time.perf_counter()-start,'passed':True}
(root/'results'/'exact_verification.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
