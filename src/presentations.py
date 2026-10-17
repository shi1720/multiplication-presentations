"""Nearest-point certificates for monomial-curve multiplication presentations.

Indices of ambient coordinates are 1..d; numpy indices are 0..d-1.
"""
from itertools import product
import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import brentq
import sympy as sp


def splits(d, mode='chain'):
    if d < 2:
        raise ValueError('d must be at least two')
    return tuple((1,k-1) if mode=='chain' else (k//2,k-k//2)
                 for k in range(2,d+1))


def all_presentations(d):
    return product(*(tuple((i,k-i) for i in range(1,k//2+1))
                     for k in range(2,d+1)))


def certificate_matrix(u, presentation):
    """H = I - 2 sum_{k=2}^d u_k A_{i_k,j_k}. Requires u_1=0."""
    u=np.asarray(u,dtype=float)
    if len(presentation)!=len(u)-1 or abs(u[0])>1e-14:
        raise ValueError('normal data and d-1 splits required')
    H=np.eye(len(u))
    for k,(i,j) in enumerate(presentation,start=2):
        if i+j!=k or not (1<=i<=j<k):
            raise ValueError('invalid split')
        if i==j: H[i-1,i-1]-=2*u[k-1]
        else:
            H[i-1,j-1]-=u[k-1]
            H[j-1,i-1]-=u[k-1]
    return H


def exact_matrix(u,presentation):
    H=sp.eye(len(u))
    for k,(i,j) in enumerate(presentation,start=2):
        if i==j: H[i-1,i-1]-=2*u[k-1]
        else:
            H[i-1,j-1]-=u[k-1]; H[j-1,i-1]-=u[k-1]
    return H


def rational_psd(M):
    """Exact Schur elimination, including singular zero-pivot cases."""
    M=sp.Matrix(M)
    while M.rows:
        p=M[0,0]
        if p<0: return False
        if p==0:
            if any(M[0,j]!=0 for j in range(1,M.cols)): return False
            M=M[1:,1:]
        else: M=M[1:,1:]-M[1:,0]*M[0,1:]/p
    return True


def face_dimension(presentation):
    return sum(i>1 for i,j in presentation)


def distance_polynomial(u):
    d=len(u); c=np.zeros(2*d+1)
    for j in range(1,d+1): c[2*j]+=1; c[j]-=2*float(u[j-1])
    return Polynomial(c)


def normal_boundary(tail):
    """Floating global minimum over all real stationary roots, not a proof oracle.
    tail=(u_3,...,u_d); returns largest allowed u_2 and minimizing parameter.
    """
    d=len(tail)+2; c=np.zeros(2*d-1); c[0]=1
    for j in range(2,d+1): c[2*j-2]+=1
    for k,uk in enumerate(tail,start=3): c[k-2]-=2*uk
    p=Polynomial(c)
    rr=p.deriv().roots()
    real=[r.real for r in rr if abs(r.imag)<1e-7]
    if not real: raise ArithmeticError('no real stationary root')
    t=min(real,key=p)
    return float(p(t)/2),float(t)


def w_boundary(d,w):
    """Stable monotone scalar solve for u3=w and remaining tail zero."""
    if w==0: return .5,0.
    sign=np.sign(w); w=abs(w)
    def deriv(s): return sum(r*s**(2*r-1) for r in range(1,d))-w
    hi=1.
    while deriv(hi)<0: hi*=2
    s=brentq(deriv,0,hi,xtol=1e-14)
    b=(1-sum((2*r-1)*s**(2*r) for r in range(1,d)))/2
    return b,float(sign*s)


def exact_nonnegative(u):
    """Exact rational real-root parity test on distance polynomial / t^2.
    Only normal data u1=0 are accepted. Positive leading coefficient is fixed.
    """
    if u[0]!=0: raise ValueError('u1 must be zero')
    t=sp.Symbol('t'); d=len(u)
    p=sp.Poly(sum(t**(2*j-2) for j in range(1,d+1))-
              2*sum(u[k-1]*t**(k-2) for k in range(2,d+1)),t)
    # A positive-leading real polynomial is nonnegative iff all real roots
    # have even multiplicity. count_roots uses exact Sturm methods over QQ.
    for f,m in sp.factor_list(p)[1]:
        if m%2 and f.count_roots(-sp.oo,sp.oo)>0: return False
    return True


def closure_boundary(tail):
    """Full quadratic closure in a d x d centered Gram matrix; CVXPY optional."""
    import cvxpy as cp
    d=len(tail)+2; H=cp.Variable((d,d),symmetric=True); v=cp.Variable()
    constraints=[H>>0]
    for s in range(2,2*d+1):
        expr=sum(H[i-1,j-1] for i in range(1,d+1)
                 for j in range(1,d+1) if i+j==s)
        rhs=(1. if s%2==0 else 0.)
        if s==2: rhs-=2*v
        elif 3<=s<=d: rhs-=2*float(tail[s-3])
        constraints.append(expr==rhs)
    problem=cp.Problem(cp.Maximize(v),constraints)
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,
                  tol_gap_rel=1e-9,max_iter=300)
    if problem.status not in ['optimal','optimal_inaccurate']:
        raise ArithmeticError(problem.status)
    hv=H.value
    residual=max(abs(sum(hv[i-1,j-1] for i in range(1,d+1)
                         for j in range(1,d+1) if i+j==s)-
                     ((1 if s%2==0 else 0)-
                     (2*v.value if s==2 else (2*tail[s-3] if 3<=s<=d else 0))))
                 for s in range(2,2*d+1))
    return {'boundary':float(v.value),'status':problem.status,
            'min_eigenvalue':float(np.linalg.eigvalsh(hv)[0]),
            'coefficient_residual':float(residual)}
