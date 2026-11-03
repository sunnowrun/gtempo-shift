"""Independent operator tests and leading-error coefficients, CPU server only."""
from foundation import *


def source_generator(hb,fs):
    p={0:hb.copy()}
    for j,f in enumerate(fs):p[1<<(2*j+1)]=f;p[1<<(2*j)]=-f.conj().T
    return p


def system(r):
    ds=1<<r;ns=[core.annihilator(r,j).conj().T@core.annihilator(r,j) for j in range(r)]
    return -.31*ns[0]+(-.11*ns[1]+1.1*ns[0]@ns[1] if r==2 else 0)


def local_audit():
    rows=[]
    for r,shared in ((1,False),(2,False),(2,True)):
        couplings=[[.7+.13j,-.22]] if r==1 else ([[.7+.13j,0],[0,.41-.09j]] if not shared else [[.7+.13j,-.22],[.19j,.41-.09j]])
        hb,fs,B=core.model([.43,-.61],couplings);hs=system(r);A=np.kron(hs,np.eye(len(hb)));H=A+B
        L=source_generator(hb,fs);L2=core.multiply(L,L,bath=True)
        Gq=np.kron(np.eye(1<<r),sum(f@f.conj().T for f in fs))
        Gn=np.kron(np.eye(1<<r),sum(f.conj().T@f for f in fs))
        for s in (0.,.5,1.):
            defect=reconstruct(L2,r,s)-B@B+(1-s)*Gq+s*Gn
            assert norm(defect)<1e-12
            rows.append(dict(r=r,shared=shared,test='square',s=s,h=0.,error=norm(defect),predicted_error=0.))
        As=normal_symbol(hs,r);poly=core.multiply({m:z*np.eye(len(hb)) for m,z in As.items()},L,bath=True)
        C=(A@B+B@A)/2-reconstruct(poly,r,1)
        Dlink=Gn/2+C
        for h in (.08,.04,.02,.01):
            U=expm(-1j*h*H)
            for method,R,pred in (
                ('q',reconstruct(core.source_cell(1j*h,hb,fs,False),r),Gq/2),
                ('normal',reconstruct(core.source_cell(1j*h,hb,fs,False),r,1),Gn/2),
                ('link',link_step(1j*h,hs,hb,fs),Dlink),
                ('weyl',q_step(1j*h,hs,hb,fs,.5,True,True),np.zeros_like(H))):
                target=expm(-1j*h*B) if method in ('q','normal') else U
                e=norm(R-target);res=norm((R-target)/h**2-pred)
                rows.append(dict(r=r,shared=shared,test=method,s=-1.,h=h,error=e,predicted_error=res))
        # Adjoint compatibility is NOT assumed for normal LINK reconstruction.
        rows.append(dict(r=r,shared=shared,test='link_adjoint',s=-1.,h=.02,
            error=norm(link_step(-.02j,hs,hb,fs)-link_step(.02j,hs,hb,fs).conj().T),predicted_error=0.))
    writecsv('local_algebra.csv',rows)


def raw_operator(lam,hs,hb,fs,N,endpoints=True):
    r=len(fs);R=reconstruct(core.source_cell(lam,hb,fs,False),r)
    A=np.kron(expm(-lam*hs),np.eye(len(hb)))
    return np.linalg.matrix_power(R@A,N)@(R if endpoints else np.eye(len(R)))


def density_obs(P,rho,O):
    X=P@rho@P.conj().T;return np.trace(O@X)/np.trace(X),np.trace(X)


def leading_error():
    rows=[]
    for name,ed,om,v,T,n0,nb in (
        ('zero',0.,0.,.7,1.,1.,0.),
        ('main',.7,1.,np.sqrt(.5),.4,1/(1+np.exp(.7)),1/(1+np.exp(1.))),
        ('filled',.7,-1.,np.sqrt(.5),.4,1/(1+np.exp(.7)),1/(1+np.exp(-1.)))):
        hb,fs,B=core.model([om],[[v]]);hs=np.diag([0.,ed]);A=np.kron(hs,np.eye(2));H=A+B
        rho=np.kron(np.diag([1-n0,n0]),np.diag([1-nb,nb]));d=np.kron(core.annihilator(1,0),np.eye(2));O=d.conj().T@d
        U=expm(-1j*T*H);Gq=np.kron(np.eye(2),fs[0]@fs[0].conj().T)
        Gn=np.kron(np.eye(2),fs[0].conj().T@fs[0]);X=density_obs(U,rho,O)[0]
        As=normal_symbol(hs,1);L=source_generator(hb,fs)
        M=reconstruct(core.multiply({m:z*np.eye(2) for m,z in As.items()},L,bath=True),1,1)
        Dlink=Gn/2+(A@B+B@A)/2-M
        Qbulk=expm_frechet(-1j*T*H,T*Gq/2,compute_expm=False)
        Qend=-.5j*(B@U+U@B)
        Nbulk=expm_frechet(-1j*T*H,T*Dlink,compute_expm=False)
        def normalized_coefficient(E):
            R=E@rho@U.conj().T+U@rho@E.conj().T
            return (np.trace(O@R)-X*np.trace(R)).real
        cq,ce,cn=map(normalized_coefficient,(Qbulk,Qend,Nbulk))
        for N in (16,32,64,128,256,512):
            h=T/N
            Rhalf=reconstruct(core.source_cell(.5j*h,hb,fs,False),1)
            Rfull=reconstruct(core.source_cell(1j*h,hb,fs,False),1)
            UA=np.kron(expm(-1j*h*hs),np.eye(2))
            Ps={'raw':raw_operator(1j*h,hs,hb,fs,N),
                'raw_noextra':np.linalg.matrix_power(q_step(1j*h,hs,hb,fs,symmetric=True),N),
                'half_endpoints':Rhalf@UA@np.linalg.matrix_power(Rfull@UA,N-1)@Rhalf,
                'link':np.linalg.matrix_power(link_step(1j*h,hs,hb,fs),N),
                'weyl':np.linalg.matrix_power(q_step(1j*h,hs,hb,fs,.5,True,True),N)}
            for method,P in Ps.items():
                obs,Z=density_obs(P,rho,O);pred=cq+ce if method=='raw' else cq if method in ('raw_noextra','half_endpoints') else cn if method=='link' else 0.
                rows.append(dict(model=name,method=method,N=N,h=h,error=(obs-X).real,error_over_h=((obs-X)/h).real,
                    predicted=pred,raw_bulk=cq,raw_endpoint=ce,link_bulk=cn,Z=Z.real,operator_error=norm(P-U)))
    writecsv('leading_error.csv',rows)


def continuum_convergence():
    rows=[]
    for r in (1,2):
        hb,fs,B=core.model([.43,-.61],[[.7+.13j,-.22]] if r==1 else [[.7+.13j,0],[0,.41-.09j]])
        hs=system(r);db=len(hb);A=np.kron(hs,np.eye(db));H=A+B;T=.8;beta=.9
        d=np.kron(core.annihilator(r,0),np.eye(db));O=d.conj().T@d
        ps=[.65,.35] if r==1 else [.4,.1,.2,.3]
        rhob=expm(-beta*hb);rhob/=np.trace(rhob);rho=np.kron(np.diag(ps),rhob)
        U=expm(-1j*T*H);thermal=expm(-beta*H);Z=np.trace(thermal)
        refs={'real':density_obs(U,rho,O)[0], 'imag':np.trace(O@thermal)/Z,
              'mixed':np.trace(U.conj().T@d@U@d.conj().T@thermal)/Z}
        for N in (8,16,32,64,128):
            for method in ('q_noextra','link','weyl'):
                def step(lam):
                    return link_step(lam,hs,hb,fs) if method=='link' else q_step(lam,hs,hb,fs,.5 if method=='weyl' else 0.,True,method=='weyl')
                P=np.linalg.matrix_power(step(1j*T/N),N);W=np.linalg.matrix_power(step(beta/N),N)
                vals={'real':density_obs(P,rho,O)[0], 'imag':np.trace(O@W)/np.trace(W),
                      'mixed':np.trace(P.conj().T@d@P@d.conj().T@W)/np.trace(P@W@P.conj().T)}
                for axis,value in vals.items():rows.append(dict(r=r,method=method,axis=axis,N=N,error=abs(value-refs[axis]),value_re=value.real,value_im=value.imag))
    writecsv('operator_convergence.csv',rows)


def green_leading():
    rows=[]
    for tag,ed,om,v,beta,T in (('main',.7,1.,np.sqrt(.5),1.,.4),('zero',0.,0.,.7,1.,1.)):
        hb,fs,B=core.model([om],[[v]]);hs=np.diag([0.,ed]);A=np.kron(hs,np.eye(2));H=A+B
        d=np.kron(core.annihilator(1,0),np.eye(2));dd=d.conj().T
        rhoS=np.diag([1.,0.]) if tag=='zero' else np.diag([1/(1+np.exp(-ed)),1/(1+np.exp(ed))])
        rhoB=np.diag([1.,0.]) if tag=='zero' else expm(-beta*hb)/np.trace(expm(-beta*hb))
        rho=np.kron(rhoS,rhoB);Gq=np.kron(np.eye(2),fs[0]@fs[0].conj().T);Gn=np.kron(np.eye(2),fs[0].conj().T@fs[0])
        L=source_generator(hb,fs);As=normal_symbol(hs,1)
        Dn=Gn/2+(A@B+B@A)/2-reconstruct(core.multiply({m:z*np.eye(2) for m,z in As.items()},L,bath=True),1,1)
        U=expm(-1j*T*H)
        for N in (16,32,64,128):
            h=T/N;R=reconstruct(core.source_cell(1j*h,hb,fs,False),1);W=R@np.kron(expm(-1j*h*hs),np.eye(2));P=np.linalg.matrix_power(W,N)@R
            V=link_step(1j*h,hs,hb,fs);Pn=np.linalg.matrix_power(V,N)
            Z=np.trace(P@rho@P.conj().T);Zn=np.trace(Pn@rho@Pn.conj().T)
            for k in range(N+1):
                t=k*h;Ut=expm(-1j*t*H);Uf=expm(-1j*(T-t)*H)
                g=np.trace(Uf@d@Ut@dd@rho@U.conj().T)
                def coefficient(Et,Ef,ET):
                    dz=np.trace(ET@rho@U.conj().T+U@rho@ET.conj().T)
                    return np.trace((Ef@d@Ut+Uf@d@Et)@dd@rho@U.conj().T+Uf@d@Ut@dd@rho@ET.conj().T)-g*dz
                def bulk(D):
                    return coefficient(expm_frechet(-1j*t*H,t*D,compute_expm=False),expm_frechet(-1j*(T-t)*H,(T-t)*D,compute_expm=False),expm_frechet(-1j*T*H,T*D,compute_expm=False))
                cq=bulk(Gq/2);cn=bulk(Dn)
                ce=coefficient(-.5j*(B@Ut+Ut@B),-.5j*(B@Uf-Uf@B),-.5j*(B@U+U@B))
                ga=np.trace(np.linalg.matrix_power(W,N-k)@d@np.linalg.matrix_power(W,k)@R@dd@rho@P.conj().T)/Z
                gn=np.trace(np.linalg.matrix_power(V,N-k)@d@np.linalg.matrix_power(V,k)@dd@rho@Pn.conj().T)/Zn
                rows.append(dict(model=tag,N=N,k=k,t=t,reference_re=g.real,reference_im=g.imag,raw_re=ga.real,raw_im=ga.imag,link_re=gn.real,link_im=gn.imag,
                    qbulk_re=cq.real,qbulk_im=cq.imag,endpoint_re=ce.real,endpoint_im=ce.imag,linkcoef_re=cn.real,linkcoef_im=cn.imag,
                    raw_coefficient_residual=abs((ga-g)/h-cq-ce),link_coefficient_residual=abs((gn-g)/h-cn)))
    writecsv('green_leading.csv',rows)


if __name__=='__main__':
    local_audit();leading_error();continuum_convergence();green_leading();print('operator audits completed')
