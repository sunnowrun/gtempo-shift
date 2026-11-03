"""Generate explicit same-target data for current Julia GMPS verification."""
from foundation import *


def jl(z):
    if isinstance(z,(bool,np.bool_)):return 'true' if z else 'false'
    if isinstance(z,(int,np.integer)):return str(z)
    z=complex(z);return f'({z.real:.17g}+({z.imag:.17g})im)'


def jmat(A):return '['+';'.join(' '.join(jl(z) for z in row) for row in A)+']'


def build_case(axis,N,method,r=1,green=False):
    if r==1:
        ed=.7;beta=1.;T=.4 if axis=='real' else .2
        hb,fs,B=core.model([1.],[[np.sqrt(.5)]]);hs=np.diag([0.,ed])
        rhoS=np.diag([1/(1+np.exp(-ed)),1/(1+np.exp(ed))])
    else:
        beta=.9;T=.8;hb,fs,B=core.model([.43,-.61],[[.7+.13j,0],[0,.41-.09j]])
        n1=core.annihilator(2,0).conj().T@core.annihilator(2,0);n2=core.annihilator(2,1).conj().T@core.annihilator(2,1)
        hs=-.31*n1-.11*n2+1.1*n1@n2;rhoS=np.diag([.4,.1,.2,.3])
    ds=len(hs);db=len(hb);H=np.kron(hs,np.eye(db))+B
    d=core.annihilator(r,0);O=d.conj().T@d;Dop=np.kron(d,np.eye(db));Op=np.kron(O,np.eye(db))
    if axis=='imag':lams=[beta/N]*N;insert={N:O};W=np.eye(db,dtype=complex);rho=None
    elif axis=='real':lams=[1j*T/N]*N+[-1j*T/N]*N;insert={N:O};W=expm(-beta*hb);W/=np.trace(W);rho=rhoS
    else:
        lams=[beta/N]*N+[1j*T/N]*N+[-1j*T/N]*N;insert={2*N:O};W=np.eye(db,dtype=complex);rho=None
        if green:insert={N:d.conj().T,2*N:d}
    cells,links=core.sequence(lams,[hs]*len(lams),hb,fs,insert,rho,False,True)
    _,links0=core.sequence(lams,[hs]*len(lams),hb,fs,{},rho,False,True)
    scalar,delta=core.quadratic(cells,W,r)
    def direct(links):
        out=np.kron(links[0],W)
        for j,c in enumerate(cells):out=np.kron(links[j+1],np.eye(db))@reconstruct(c,r,.5 if method=='weyl' else 0.)@out
        if method=='weyl':out*=np.prod([weyl.coupling_scalar(z,fs) for z in lams])
        return np.trace(out)
    num,den=direct(links),direct(links0)
    U=expm(-1j*T*H);thermal=expm(-beta*H)
    if axis=='real':physical=density_value(U,np.kron(rhoS,W),Op)
    elif axis=='imag':physical=np.trace(Op@thermal)/np.trace(thermal)
    else:physical=np.trace(U.conj().T@(Dop if green else Op)@U@(Dop.conj().T if green else np.eye(len(H)))@thermal)/np.trace(thermal)
    scale=np.prod([weyl.coupling_scalar(z,fs) for z in lams]) if method=='weyl' else 1.
    # Full polynomial is an independent oracle only for the small cases.
    K=core.K_polynomial(links,r) if N<=2 else {}
    K0=core.K_polynomial(links0,r) if r==2 else {}
    if method=='weyl':
        for j in range(r,(len(cells)+1)*r):K=contact(K,j)
        for j in range(r,(len(cells)+1)*r):K0=contact(K0,j)
    if K:
        IF=core.exponential(core.quadratic_poly(delta),len(cells)*r)
        IF={m<<(2*r):scalar*z for m,z in IF.items()}
        P=core.multiply(K,IF)
        for j in range((len(cells)+1)*r):P=core.multiply(P,{0:1,3<<(2*j):1})
        val=scale*P.get((1<<(2*(len(cells)+1)*r))-1,0)
        assert abs(val-num)<1e-10,(axis,N,r,method,val,num)
    return dict(name=f'{axis}_{method}_r{r}_N{N}_g{int(green)}',axis=axis,N=N,r=r,L=len(cells),green=green,method=method,
        delta=delta,links=links,links0=links0,scalar=scalar,scale=scale,num=num,den=den,physical=physical,K=K,K0=K0)


def density_value(U,rho,O):
    X=U@rho@U.conj().T;return np.trace(O@X)/np.trace(X)


def contact(K,pair):
    b=3<<(2*pair);out={}
    for m,z in K.items():
        k=(m>>(2*pair))&3;base=m&~b
        if k==0:
            out[base]=out.get(base,0)+.5*z;out[base|b]=out.get(base|b,0)+.5*z
        elif k==3:
            out[base]=out.get(base,0)-.5*z;out[base|b]=out.get(base|b,0)+1.5*z
        else:out[m]=out.get(m,0)+z
    return out


def run():
    cases=[]
    for axis in ('imag','real','mixed'):
        for method in ('q','weyl'):
            for N in (2,4,8,16):cases.append(build_case(axis,N,method))
    for axis in ('imag','real','mixed'):
        for method in ('q','weyl'):cases.append(build_case(axis,1,method,r=2))
    for N in (2,4,8,16):cases.append(build_case('mixed',N,'weyl',green=True))
    with open(ROOT/'results/package_cases.jl','w') as f:
        f.write('const CASES = [\n')
        for c in cases:
            fields=[f'name="{c["name"]}"',f'axis=:{c["axis"]}',f'method=:{c["method"]}']
            fields += [f'{k}={jl(c[k])}' for k in ('N','r','L','green','scalar','scale','num','den','physical')]
            fields += ['delta='+jmat(c['delta'])]
            for k in ('links','links0'):fields.append(k+'=['+','.join(jmat(A) for A in c[k])+']')
            fields.append('kmasks=UInt64['+','.join(str(m) for m in c['K'])+']')
            fields.append('kcoeff=ComplexF64['+','.join(jl(z) for z in c['K'].values())+']')
            fields.append('k0masks=UInt64['+','.join(str(m) for m in c['K0'])+']')
            fields.append('k0coeff=ComplexF64['+','.join(jl(z) for z in c['K0'].values())+']')
            f.write('(;'+','.join(fields)+'),\n')
        f.write(']\n')
    print('generated',len(cases),'cases; small full Grassmann integrals passed')


if __name__=='__main__':run()
