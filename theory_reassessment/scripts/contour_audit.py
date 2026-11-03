"""Derive raw and link IF coefficients from finite source exponentials."""
from foundation import *
from operator_audit import raw_operator,source_generator


def source_if(case,nodes,links=False,half=False):
    axis=case['axis'];N=int(case['N']);Ni=int(case['Ni']);T=float(case['T']);beta=float(case['beta'])
    hb,fs,B=core.model([float(case['omega'])],[[np.sqrt(float(case['weight']))]])
    seq=[];row=[];col=[];sg=[]
    if axis!='real':
        for j in range(1,Ni+1):
            seq.append(beta/Ni);row.append(('τ',j+1));col.append(('τ',j if links else j+1));sg.append(-1)
    if axis!='imag':
        count=N if links else N+1
        for j in range(1,count+1):
            seq.append(1j*T/N*(.5 if half and j in (1,count) else 1));row.append(('+',j+1 if links else j));col.append(('+',j));sg.append(1)
        for j in range(count,0,-1):
            seq.append(-1j*T/N*(.5 if half and j in (1,count) else 1));row.append(('-',j));col.append(('-',j+1 if links else j));sg.append(-1 if axis=='real' else 1)
    cache={};cells=[]
    for z,s in zip(seq,sg):
        if (z,s) not in cache:
            c=core.source_cell(z,hb,fs,False)
            cache[z,s]={m:v*s**m.bit_count() for m,v in c.items()}
        cells.append(cache[z,s])
    W=expm(-beta*hb) if axis=='real' else np.eye(len(hb),dtype=complex)
    W/=np.trace(expm(-beta*hb))
    weight,D=core.quadratic(cells,W,1)
    # A transpose at backward thermal trace is never assumed: all signs derived.
    target=np.zeros((len(nodes),len(nodes)),complex);ix={p:i for i,p in enumerate(nodes)}
    for i,r in enumerate(row):
        for j,c in enumerate(col):target[ix[r],ix[c]]+=D[i,j]
    return target,weight


def compare_sources():
    rows=[];details=[]
    for case in csv.DictReader(open(ROOT/'results/cases.csv')):
        if not case['id'].startswith(('main','zero','filled')) or int(case['N']) not in (2,4,8,16):continue
        nodes,ds=load_case(case);D,weight=source_if(case,nodes)
        # Record per block to expose signs/contact mismatches, not only norms.
        for b in sorted({b for b,j in nodes}):
            for c in sorted({b for b,j in nodes}):
                ix=[i for i,p in enumerate(nodes) if p[0]==b];jx=[i for i,p in enumerate(nodes) if p[0]==c]
                A=D[np.ix_(ix,jx)];B=ds['raw'][np.ix_(ix,jx)]
                rows.append(dict(id=case['id'],b1=b,b2=c,error=norm(A-B),plus_error=norm(A+B),norm=norm(A),weight_re=weight.real))
        if int(case['N'])==4:
            np.savez(ROOT/'results'/f"source_{case['id']}.npz",raw=D,upstream=ds['raw'],link=source_if(case,nodes,True)[0],vertex=ds['vertex'])
        g,Z,C=gaussian(case,nodes,D)
        ref=physical_reference(case)
        details.append(dict(id=case['id'],error=norm(g-ref)/norm(ref),g0=g[0].real,Z=Z.real))
    writecsv('source_blocks.csv',rows);writecsv('source_physics.csv',details)


def wrong_continuum(case):
    ed,om,v2,beta,T=[float(case[k]) for k in ('ed','omega','weight','beta','T')]
    hb,fs,B=core.model([om],[[np.sqrt(v2)]]);hs=np.diag([0.,ed]);A=np.kron(hs,np.eye(2));H=A+B
    Hflip=A+np.kron(np.eye(2),hb)-(B-np.kron(np.eye(2),hb))
    R=expm(-beta*H);d=np.kron(core.annihilator(1,0),np.eye(2));N=int(case['N'])
    # The two real branches have different Hamiltonians, so future evolution
    # does NOT cancel for an internal insertion.
    U=expm(-1j*T*Hflip);V=expm(1j*T*H);Z=np.trace(V@U@R)
    return np.array([np.trace(V@expm(-1j*(T-t)*Hflip)@d@expm(-1j*t*Hflip)@d.conj().T@R)/Z for t in np.linspace(0,T,N+1)]),Z/np.trace(expm(-beta*hb))


def link_and_sign_audit():
    rows=[]
    for case in csv.DictReader(open(ROOT/'results/cases.csv')):
        if not case['id'].startswith(('main','short')):continue
        axis=case['axis'];N=int(case['N'])
        nodes,ds=load_case(case)
        # Quadratic source extraction is O(N²); no exponential exterior object.
        link,weight=source_if(case,nodes,True)
        if axis=='real':
            e=norm(link-ds['vertex']);assert e<1e-10
        else:e=norm(link-ds['vertex'])
        ref=physical_reference(case)
        variants={'derived_link':link}
        if axis!='imag':variants['half_endpoints']=source_if(case,nodes,half=True)[0]
        if axis=='mixed':
            no=ds['vertex'].copy()
            for i,(b,j) in enumerate(nodes):
                for k,(c,l) in enumerate(nodes):
                    no[i,k]*=(-1)**((b=='+' and j>=2)^(c=='+'))
            variants['vertex_nosign']=no
            variants['raw_flip_plus']=ds['raw']*np.array([-1 if b=='+' else 1 for b,j in nodes])[:,None]*np.array([-1 if b=='+' else 1 for b,j in nodes])[None,:]
            wrong,Zwrong=wrong_continuum(case)
            g,Z,C=gaussian(case,nodes,ds['vertex'])
            rows.append(dict(id=case['id'],axis=axis,N=N,method='old_vertex_vs_wrong_limit',error=norm(g-wrong)/norm(wrong),physical_error=norm(g-ref)/norm(ref),kernel_error=e,Z_error=abs(Z-Zwrong)))
        for method,D in variants.items():
            g,Z,C=gaussian(case,nodes,D)
            rows.append(dict(id=case['id'],axis=axis,N=N,method=method,error=norm(g-ref)/norm(ref),physical_error=norm(g-ref)/norm(ref),kernel_error=e,Z_error=0.))
    writecsv('link_and_sign.csv',rows)


if __name__=='__main__':compare_sources();link_and_sign_audit()
