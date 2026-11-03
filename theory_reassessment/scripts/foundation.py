"""Finite-bath reconstruction and independent Gaussian contour integration.

All numerical execution belongs on the CPU server. External reference utilities
were copied verbatim; their provenance is recorded by source_audit.py.
"""
from pathlib import Path
import sys, csv, json
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
import core
import weyl
import numpy as np
from scipy.linalg import expm, expm_frechet, norm


def reconstruct(poly,r,s=0.):
    """s=0: original Q; s=1/2: Weyl; s=1: normal coherent LINK symbol."""
    ds=1<<r; db=len(poly[0]); out=np.zeros((ds*db,ds*db),complex)
    I=np.eye(ds);P=core.parity(r)
    for mask,c in poly.items():
        q=I.copy()
        for j in range(r):
            d=core.annihilator(r,j);k=mask>>(2*j)&3
            q=q@(I if k==0 else d if k==1 else d.conj().T if k==2 else (1-s)*I-d.conj().T@d)
        out+=np.kron(q@(P if mask.bit_count()%2 else I),c)
    return out


def normal_symbol(G,r):
    # Audited for arbitrary one mode and diagonal interacting multi-mode systems.
    k=core.kernel(G,0,0,r)
    overlap={}
    for j in range(r):overlap.update(core.word([2*j+1,2*j],-1))
    return core.multiply(core.exponential(overlap,r),k)


def link_step(lam,hs,hb,fs):
    r=len(fs);cell=core.source_cell(lam,hb,fs,False)
    sym=normal_symbol(expm(-lam*hs),r)
    sym={m:z*np.eye(len(hb)) for m,z in sym.items()}
    return reconstruct(core.multiply(sym,cell,bath=True),r,1.)


def q_step(lam,hs,hb,fs,s=0.,symmetric=False,compensate=False):
    r=len(fs);db=len(hb);R=reconstruct(core.source_cell(lam,hb,fs,False),r,s)
    if compensate:
        c=sum(np.trace(f@f.conj().T+f.conj().T@f)/db for f in fs)
        R*=np.exp(lam*lam*c/4)
    A=np.kron(expm(-lam*hs/(2 if symmetric else 1)),np.eye(db))
    return A@R@A if symmetric else R@A


def kernel_matrix(case,nodes):
    axis=case['axis'];N=int(case['N']);Ni=int(case['Ni'])
    dt=float(case['T'])/N;ht=float(case['beta'])/Ni;ed=float(case['ed'])
    n=len(nodes);K=np.zeros((n,n),complex);ix={p:i for i,p in enumerate(nodes)}
    def put(out,inn,z):K[ix[out],ix[inn]]+=z
    if axis!='real':
        for j in range(1,Ni+1):put(('τ',j+1),('τ',j),np.exp(-ht*ed))
    if axis!='imag':
        for j in range(1,N+1):
            put(('+',j+1),('+',j),np.exp(-1j*dt*ed))
            put(('-',j),('-',j+1),np.exp(1j*dt*ed))
    scalar=1.
    if axis=='imag':put(('τ',1),('τ',Ni+1),-1)
    elif axis=='real':
        put(('-',N+1),('+',N+1),-1)
        a=np.exp(-float(case['beta'])*ed)
        put(('+',1),('-',1),a);scalar=1/(1+a)
    else:
        put(('-',N+1),('+',N+1),1)
        put(('τ',1),('-',1),1)
        put(('+',1),('τ',Ni+1),-1)
    return K,scalar


def load_case(case):
    id=case['id'];nodes=[]
    for b,j in csv.reader(open(ROOT/'results'/f'{id}_nodes.csv')):nodes.append((b,int(j)))
    ds={}
    for mode in ('raw','vertex'):
        arr=np.loadtxt(ROOT/'results'/f'{id}_{mode}.csv',delimiter=',',ndmin=2)
        ds[mode]=arr[:,:len(nodes)]+1j*arr[:,len(nodes):]
    return nodes,ds


def gaussian(case,nodes,D):
    K,scalar=kernel_matrix(case,nodes);M=K+D-np.eye(len(K));C=-np.linalg.inv(M)
    Z=np.linalg.det(-M)*scalar
    axis=case['axis'];b='τ' if axis=='imag' else '+';ix={p:i for i,p in enumerate(nodes)}
    ids=[i for i,(br,j) in enumerate(nodes) if br==b]
    g=C[ids,ix[b,1]]
    return g,Z,C


def physical_reference(case):
    ed,om,v2,beta,T=[float(case[k]) for k in ('ed','omega','weight','beta','T')]
    N=int(case['N']);Ni=int(case['Ni']);axis=case['axis'];H=np.array([[ed,np.sqrt(v2)],[np.sqrt(v2),om]])
    f=np.linalg.inv(np.eye(2)+expm(beta*H));n0=1/(1+np.exp(beta*ed))
    if axis=='imag':return np.array([(expm(-t*H)@(np.eye(2)-f))[0,0] for t in np.linspace(0,beta,Ni+1)])
    return np.array([(expm(-1j*t*H)@(np.eye(2)-f) if axis=='mixed' else expm(-1j*t*H)*(1-n0))[0,0] for t in np.linspace(0,T,N+1)])


def writecsv(name,rows):
    with open(ROOT/'results'/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def gaussian_audit():
    rows=[]
    for case in csv.DictReader(open(ROOT/'results/cases.csv')):
        nodes,ds=load_case(case);ref=physical_reference(case)
        for mode,D in ds.items():
            g,Z,C=gaussian(case,nodes,D)
            rows.append(dict(id=case['id'],axis=case['axis'],N=case['N'],Ni=case['Ni'],method=mode,
                error=norm(g-ref)/norm(ref),g0_re=g[0].real,g0_im=g[0].imag,
                gend_re=g[-1].real,gend_im=g[-1].imag,Z_re=Z.real,Z_im=Z.imag))
    writecsv('gaussian.csv',rows)


if __name__=='__main__':gaussian_audit()
