"""Weyl GTEMPO, finite-bath reference and explicit repeated local tensors.
Parent code is imported read-only as an independent exterior-algebra utility.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import core as old
import numpy as np
from scipy.linalg import expm,norm


def gmul(p,q,left_parity=0):
    out={}
    for a,x in p.items():
        for b,y in q.items():
            if a&b:continue
            s=old.sign(a,b)*(-1)**((a.bit_count()+left_parity)*b.bit_count())
            out[a|b]=out.get(a|b,0)+s*(x@y)
    return out


def weyl_kernel(r,gaussian_measure=False):
    ds=1<<r;I=np.eye(ds,dtype=complex);p={0:I}
    for j in range(r):
        d=old.annihilator(r,j);n=d.conj().T@d
        local={0:I/2-n,1<<(2*j):-d.conj().T,1<<(2*j+1):d,
               3<<(2*j):(I/2+n if gaussian_measure else I)}
        p=gmul(p,local)
    return p


def weyl_reconstruct(poly,r):
    ds=1<<r;db=len(poly[0]);out=np.zeros((ds*db,ds*db),complex)
    I=np.eye(ds)
    for mask,c in poly.items():
        q=I.copy()
        for j in range(r):
            d=old.annihilator(r,j);k=(mask>>(2*j))&3
            q=q@(I if k==0 else d if k==1 else d.conj().T if k==2 else I/2-d.conj().T@d)
        out+=np.kron(q@(old.parity(r) if mask.bit_count()%2 else I),c)
    return out


def coupling_scalar(lam,fs):
    c=sum(np.trace(f@f.conj().T+f.conj().T@f)/len(f) for f in fs)
    return np.exp(lam*lam*c/4)


def step(lam,hs,hb,fs,method='weyl'):
    r=len(fs);db=len(hb)
    cell=old.source_cell(lam,hb,fs,corrected=method=='bath_corrected')
    R=weyl_reconstruct(cell,r)*coupling_scalar(lam,fs) if method=='weyl' else old.reconstruct(cell,r)
    A=np.kron(expm(-lam*hs/(1 if method=='order1' else 2)),np.eye(db))
    return A@R if method=='order1' else A@R@A


def setup(lams,systems,hb,fs,insertions=None,rhoS=None,method='weyl'):
    cells,links=old.sequence(lams,systems,hb,fs,insertions,rhoS,corrected=method=='bath_corrected',strang=method!='order1')
    scalar=np.prod([coupling_scalar(z,fs) for z in lams]) if method=='weyl' else 1.
    return cells,links,scalar


def dense_trace(cells,links,W,r,scalar=1.):
    db=len(W);p=np.kron(links[0],W)
    for j,c in enumerate(cells):p=np.kron(links[j+1],np.eye(db))@weyl_reconstruct(c,r)@p
    return scalar*np.trace(p)


def K_polynomial(links,r,parities=None,gaussian_measure=False):
    parities=parities or [0]*len(links);M=weyl_kernel(r,gaussian_measure)
    p={0:links[0]}
    for j in range(len(links)-1):
        p=gmul(old.shift(M,j,r),p)
        p=gmul({0:links[j+1]},p,parities[j+1])
    return {m:np.trace(x) for m,x in p.items()}


def grassmann_integral(cells,links,W,r,scalar=1.,parities=None,gaussian_measure=False):
    K=K_polynomial(links,r,parities,gaussian_measure)
    weight,D=old.quadratic(cells,W,r)
    IF={m:weight*x for m,x in old.exponential(old.quadratic_poly(D),len(cells)*r).items()}
    poly=old.multiply(K,IF)
    if gaussian_measure:
        for j in range(len(cells)*r):poly=old.multiply(poly,{0:1,3<<(2*j):1})
    return scalar*poly.get((1<<(2*len(cells)*r))-1,0)


def local_K(lam,hs,r,fs=None,gaussian_measure=False):
    A=expm(-lam*hs/2);s=1 if fs is None else coupling_scalar(lam,fs)
    return {m:s*A@x@A for m,x in weyl_kernel(r,gaussian_measure).items()}


def polynomial_from_uniform(local,N,weight=None,O=None):
    d=len(local[0]);weight=np.eye(d) if weight is None else weight;O=np.eye(d) if O is None else O
    p={0:weight}
    bits=int(np.log2(len(local)))
    for j in range(N):p=gmul({m<<(j*bits):v for m,v in local.items()},p)
    return {m:np.trace(O@x) for m,x in p.items()}


def folded_coefficients(plus,minus,N,weight,O=None):
    """Repeated folded supertensor; canonical local order (+a,+bar,-a,-bar).
    q remembers parity of minus variables from earlier physical-time cells.
    """
    r=int(np.log2(len(plus))/2);d=len(weight);O=np.eye(d) if O is None else O
    p={(0,0):weight}
    for j in range(N):
        out={}
        for (mask,q),X in p.items():
            for a,A in plus.items():
                for b,B in minus.items():
                    s=(-1)**(q*(a.bit_count()+b.bit_count()))
                    key=(mask|((a|(b<<(2*r)))<<(4*r*j)),q^(b.bit_count()%2))
                    out[key]=out.get(key,0)+s*A@X@B
        p=out
    ans={}
    for (mask,q),X in p.items():ans[mask]=ans.get(mask,0)+np.trace(O@X)
    return ans


def permute_polynomial(p,old_to_new):
    out={}
    for mask,v in p.items():
        indices=[old_to_new[i] for i in range(len(old_to_new)) if mask>>i&1]
        term=old.word(indices,v)
        for m,x in term.items():out[m]=out.get(m,0)+x
    return out


def uniform_integrated(localK,localE,N,weightS,weightB,r=1,O=None):
    """Exact repeated contraction of separate local K and IF tensors.
    Uses one parity bit for cross-cell Koszul signs; bare measure.
    """
    ds=len(weightS);db=len(weightB);D=ds*db;top=(1<<(2*r))-1
    T=np.zeros((2*D,2*D),complex)
    for p in (0,1):
        for m,A in localK.items():
            n=top^m;B=localE[n];q=p^(n.bit_count()%2)
            s=old.sign(m,n)*(-1)**(m.bit_count()*p)
            T[q*D:(q+1)*D,p*D:(p+1)*D]+=s*np.kron(A,B)
    power=np.linalg.matrix_power(T,N)[:D,:D]
    Of=np.kron(np.eye(ds) if O is None else O,np.eye(db))
    return np.trace(Of@power@np.kron(weightS,weightB))


def mixed_repeated_coefficients(imag,plus,minus,Ni,Nt,initial_insert=None,turn=None,odd_pair=False):
    """Imaginary repeated segment joined to repeated folded real segment.
    odd_pair inserts two odd system matrices without treating their coefficients
    as commuting; the resulting extra sign is (-1)^total plus degree.
    """
    r=int(np.log2(len(imag))/2);bits=2*r;d=len(imag[0]);p={0:np.eye(d)}
    for j in range(Ni):p=gmul({m<<(bits*j):x for m,x in imag.items()},p)
    X0=np.eye(d) if initial_insert is None else initial_insert
    p={(m,0,0):X0@x for m,x in p.items()}
    for j in range(Nt):
        out={}
        for (mask,q,pp),X in p.items():
            for a,A in plus.items():
                for b,B in minus.items():
                    s=(-1)**(q*(a.bit_count()+b.bit_count()))
                    key=(mask|((a|(b<<bits))<<(bits*Ni+2*bits*j)),q^(b.bit_count()%2),pp^(a.bit_count()%2))
                    out[key]=out.get(key,0)+s*A@X@B
        p=out
    O=np.eye(d) if turn is None else turn;ans={}
    for (mask,q,pp),X in p.items():ans[mask]=ans.get(mask,0)+(-1 if odd_pair and pp else 1)*np.trace(O@X)
    return ans
