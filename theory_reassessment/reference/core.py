"""Finite-bath order-two prototype. No Julia/package code is imported here."""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.linalg import expm


def annihilator(n, j):
    a=np.zeros((1<<n,1<<n),complex)
    for k in range(1<<n):
        if k>>j&1: a[k^(1<<j),k]=(-1)**((k&((1<<j)-1)).bit_count())
    return a


def parity(n): return np.diag([(-1)**k.bit_count() for k in range(1<<n)])


def sign(a,b):
    return (-1)**sum((b&((1<<j)-1)).bit_count() for j in range(a.bit_length()) if a>>j&1)


def multiply(p,q,bath=False,max_degree=None):
    out={}
    for a,x in p.items():
        for b,y in q.items():
            if a&b or (max_degree is not None and (a|b).bit_count()>max_degree): continue
            s=sign(a,b)*((-1)**(a.bit_count()*b.bit_count()) if bath else 1)
            v=s*(x@y if bath else x*y)
            out[a|b]=out.get(a|b,0)+v
    return out


def exponential(p,pairs):
    out={0:1.+0j};term=out
    for k in range(1,2*pairs+1):
        term={a:v/k for a,v in multiply(term,p).items() if abs(v)>1e-30}
        if not term: break
        for a,v in term.items(): out[a]=out.get(a,0)+v
    return out


def word(indices,value=1):
    mask=0;s=1
    for j in indices:
        bit=1<<j
        if mask&bit: return {}
        s*=sign(mask,bit);mask|=bit
    return {mask:s*value}


def add(p,q):
    out=p.copy()
    for k,v in q.items():out[k]=out.get(k,0)+v
    return out


def model(energies,couplings):
    bs=[annihilator(len(energies),j) for j in range(len(energies))]
    hb=sum(e*b.conj().T@b for e,b in zip(energies,bs))
    fs=[sum(v*b for v,b in zip(row,bs)) for row in couplings]
    r=len(fs);db=len(hb);B=np.kron(np.eye(1<<r),hb)
    for j,f in enumerate(fs):
        d=np.kron(annihilator(r,j),np.eye(db));F=np.kron(parity(r),f)
        B+=d.conj().T@F+F.conj().T@d
    return hb,fs,B


def source_cell(lam,hb,fs,corrected=True):
    """Exact finite exterior-algebra representation of exp(-lam L+lam² C/2)."""
    r=len(fs);db=len(hb);de=1<<(2*r)
    C=sum(f@f.conj().T for f in fs)
    generator=np.kron(np.eye(de),-lam*hb+(lam*lam/2*C if corrected else 0))
    pe=parity(2*r)
    for j,f in enumerate(fs):
        aa=[]
        for pos in (2*j,2*j+1):
            x=np.zeros((de,de))
            for k in range(de):
                if not (k>>pos&1):x[k|(1<<pos),k]=(-1)**((k&((1<<pos)-1)).bit_count())
            aa.append(x)
        generator-=lam*np.kron(aa[1]@pe,f)
        generator-=lam*np.kron(pe@aa[0],f.conj().T)
    e=expm(generator)
    return {a:e[db*a:db*(a+1),:db] for a in range(de)}


def reconstruct(poly,r):
    ds=1<<r;db=len(poly[0]);out=np.zeros((ds*db,ds*db),complex)
    for mask,c in poly.items():
        q=np.eye(ds,dtype=complex)
        for j in range(r):
            d=annihilator(r,j);k=mask>>(2*j)&3
            q=q@(np.eye(ds) if k==0 else d if k==1 else d.conj().T if k==2 else np.eye(ds)-d.conj().T@d)
        out+=np.kron(q@(parity(r) if mask.bit_count()%2 else np.eye(ds)),c)
    return out


def shift(poly,j,r):return {a<<(2*r*j):c for a,c in poly.items()}


def bath_polynomial(cells,W,r,max_degree=None):
    p={0:W.copy()}
    for j,cell in enumerate(cells):p=multiply(shift(cell,j,r),p,True,max_degree)
    return {a:np.trace(c) for a,c in p.items()}


def quadratic(cells,W,r):
    """O(L²) insertions; coefficients exact for the finite cell exponentials.
    Does not reconstruct the IF from only its second degree as an assumption:
    full-polynomial Wick audits are supplied separately.
    """
    L=len(cells);db=len(W)
    prefix=[np.eye(db,dtype=complex)]
    for cell in cells:prefix.append(cell[0]@prefix[-1])
    suffix=[None]*(L+1);suffix[L]=np.eye(db,dtype=complex)
    for j in range(L-1,-1,-1):suffix[j]=suffix[j+1]@cells[j][0]
    weight=np.trace(prefix[-1]@W)
    if abs(weight)<1e-13:raise ValueError('source-free Gaussian weight is singular')
    out=np.zeros((L*r,L*r),complex)
    for k in range(L):
        for alpha in range(r):
            for gamma in range(r):
                abit=1<<(2*gamma);bbit=1<<(2*alpha+1)
                c=cells[k][abit|bbit]
                out[k*r+alpha,k*r+gamma]=sign(bbit,abit)*np.trace(suffix[k+1]@c@prefix[k]@W)/weight
        middle=np.eye(db,dtype=complex)
        for l in range(k-1,-1,-1):
            for alpha in range(r):
                for gamma in range(r):
                    # bar in late cell / a in early cell
                    c=np.trace(suffix[k+1]@cells[k][1<<(2*alpha+1)]@middle@cells[l][1<<(2*gamma)]@prefix[l]@W)/weight
                    out[k*r+alpha,l*r+gamma]=-c
                    # a in late cell / bar in early cell
                    c=np.trace(suffix[k+1]@cells[k][1<<(2*gamma)]@middle@cells[l][1<<(2*alpha+1)]@prefix[l]@W)/weight
                    out[l*r+alpha,k*r+gamma]=c
            middle=middle@cells[l][0]
    return weight,out


def quadratic_poly(delta):
    p={}
    for i in range(len(delta)):
        for j in range(len(delta)):
            p=add(p,word([2*i+1,2*j],delta[i,j]))
    return p


def dense_trace(cells,links,W,r):
    db=len(W);p=np.kron(links[0],W);cache={}
    for j,c in enumerate(cells):
        if id(c) not in cache: cache[id(c)]=reconstruct(c,r)
        p=np.kron(links[j+1],np.eye(db))@cache[id(c)]@p
    return np.trace(p)


def kernel(G,outnode,innode,r,negative=False):
    """Diagonal multimode or arbitrary one-mode coherent-state kernel."""
    if r==1:
        a=2*innode;b=2*outnode+1;s=-1 if negative else 1
        p={0:G[0,0]}
        p=add(p,word([a],G[0,1]));p=add(p,word([b],s*G[1,0]));p=add(p,word([b,a],s*G[1,1]))
        return p
    assert np.linalg.norm(G-np.diag(np.diag(G)))<1e-12
    out={}
    for state,value in enumerate(np.diag(G)):
        p={0:value}
        for alpha in range(r):
            if state>>alpha&1:
                p=multiply(p,word([2*(r*outnode+alpha)+1,2*(r*innode+alpha)],-1 if negative else 1))
        out=add(out,p)
    return out


def K_polynomial(links,r):
    L=len(links)-1
    p=kernel(links[L],0,L,r,True)
    for j in range(L-1,0,-1):p=multiply(p,kernel(links[j],j+1,j,r))
    return multiply(p,kernel(links[0],1,0,r))


def grassmann_trace(cells,links,W,r,use_wick=True):
    L=len(cells);pairs=r*(L+1)
    if use_wick:
        scalar,delta=quadratic(cells,W,r)
        bath={a<<(2*r):v*scalar for a,v in exponential(quadratic_poly(delta),L*r).items()}
    else:bath={a<<(2*r):v for a,v in bath_polynomial(cells,W,r).items()}
    p=multiply(K_polynomial(links,r),bath)
    for j in range(pairs):p=multiply(p,{0:1,3<<(2*j):1})
    return p.get((1<<(2*pairs))-1,0)


def sequence(lambdas,systems,hb,fs,insertions=None,rhoS=None,corrected=True,strang=True):
    insertions=insertions or {};r=len(fs);ds=1<<r
    cells=[];halves=[];cache={}
    for lam,Hs in zip(lambdas,systems):
        key=complex(lam)
        if key not in cache:cache[key]=source_cell(lam,hb,fs,corrected)
        cells.append(cache[key]);halves.append(expm(-lam*Hs/(2 if strang else 1)))
    # chronological right-to-left factors: L_j R_j Rhalf_j;
    # uncorrected Lie mode has no right system half.
    right=halves if strang else [np.eye(ds)]*len(cells)
    links=[right[0]@(np.eye(ds) if rhoS is None else rhoS)]
    for j in range(1,len(cells)):
        links.append(right[j]@insertions.get(j,np.eye(ds))@halves[j-1])
    links.append(insertions.get(len(cells),np.eye(ds))@halves[-1])
    return cells,links


def toeplitz_defect(a):
    return float(np.max(np.abs(a[1:,1:]-a[:-1,:-1]))) if len(a)>1 else 0.
