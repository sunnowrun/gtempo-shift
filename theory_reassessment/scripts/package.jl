module TheoryPackage
using GTEMPO,QuAPI,Z2Tensors,LinearAlgebra,Statistics
using GTEMPO: @tensor
const ROOT=dirname(@__DIR__)
include(joinpath(ROOT,"results","package_cases.jl"))

# Complex-valued coefficient wrapper; no global package method is overwritten.
struct ExplicitCorr <: QuAPI.AbstractCorrelationFunction
    D::Matrix{ComplexF64}
end
QuAPI.branches(::ExplicitCorr)=(:τ,)
QuAPI.branch(c::ExplicitCorr;b1::Symbol,b2::Symbol)=c.D
function Base.getindex(c::ExplicitCorr,i::GrassmannContourIndex,j::GrassmannContourIndex)
    i.conj==j.conj && return 0.0im
    i.conj ? -c.D[i.j,j.j] : c.D[i.j,j.j]
end

function layout(c)
    # This lattice is a storage host for the explicitly ordered contour cells.
    # The physical contour, branches and steps are in the generated input.
    lat=ImagGrassmannLattice(N=c.L,bands=c.r,δτ=1.)
    positions=[(lat[GrassmannContourIndex(j,:τ,b,false)],lat[GrassmannContourIndex(j,:τ,b,true)]) for j in 1:c.L+1 for b in 1:c.r]
    lat,positions
end

function build_if(c,lat,pos,trunc)
    if c.r==1
        p=vacuumstate(ComplexF64,lat)
        return hybriddynamics!(p,lat,ExplicitCorr(c.delta);trunc)
    end
    # Existing single-band row API excludes cross bands. Here independent baths
    # imply block-diagonal flavor coefficients; build each physical row explicitly.
    p=vacuumstate(ComplexF64,lat)
    for i in axes(c.delta,1)
        ids=[j for j in axes(c.delta,2) if abs(c.delta[i,j])>1e-16]
        isempty(ids)&&continue
        row=pos[c.r+i][2];cols=[pos[c.r+j][1] for j in ids]
        f=GTEMPO.partialmpo(row,cols,c.delta[i,ids])*vacuumstate(ComplexF64,lat)
        p=mult!(p,f;trunc)[1]
    end
    p
end

function coherent_K(c,lat,pos,links,trunc)
    @assert c.r==1
    p=vacuumstate(ComplexF64,lat);factor=1.0+0im;odd=Int[]
    for j in c.L:-1:0
        out=j==c.L ? 1 : j+2;inn=j==0 ? 1 : j+1
        bar=pos[out][2];a=pos[inn][1];sg=j==c.L ? -1 : 1;G=links[j+1]
        if abs(G[1,2])+abs(G[2,1])<1e-14
            if abs(G[1,1])>1e-15
                factor*=G[1,1];apply!(exp(GTerm(bar,a;coeff=sg*G[2,2]/G[1,1])),p)
            else
                apply!(GTerm(bar,a;coeff=sg*G[2,2]),p)
            end
        elseif abs(G[1,2])>1e-14
            @assert abs(G[2,1])<1e-14
            push!(odd,a);factor*=G[1,2]
        else
            push!(odd,bar);factor*=sg*G[2,1]
        end
    end
    isempty(odd)||apply!(GTerm(Tuple(odd);coeff=1.),p)
    canonicalize!(p,alg=Orthogonalize(trunc=trunc))
    factor*p
end

function contact_pair(p,a,bar,trunc)
    x,y=minmax(a,bar);@assert y==x+1
    eta=bar>a ? 1 : -1
    np=copy(p);nq=copy(p)
    @tensor tmp[1,2;3] := GTEMPO.n̂[2,4]*np[x][1,4,3]
    np[x]=tmp;GTEMPO.unset_svectors!(np)
    @tensor tmp[1,2;3] := GTEMPO.n̂[2,4]*nq[y][1,4,3]
    nq[y]=tmp;GTEMPO.unset_svectors!(nq)
    create=convert(PartialMPO,GTerm(x,y;coeff=1.))
    remove=PartialMPO([GTEMPO.σ₋,adjoint(GTEMPO.σ₊)],[x,y])
    result=.5*p+.5*np+.5*nq+.5*eta*(create*p)-.5*eta*(remove*p)
    canonicalize!(result,alg=Orthogonalize(trunc=trunc))
    result
end

function polynomial_K(c,lat,pos,masks,coeffs,trunc)
    p=0.0*vacuumstate(ComplexF64,lat)
    for (m,z) in zip(masks,coeffs)
        abs(z)<1e-15&&continue
        inds=Int[pos[j][k+1] for j in eachindex(pos) for k in 0:1 if !iszero((m>>(2*(j-1)+k))&1)]
        f=vacuumstate(ComplexF64,lat)
        isempty(inds) ? (f=z*f) : apply!(GTerm(Tuple(inds);coeff=z),f)
        p=p+f
        canonicalize!(p,alg=Orthogonalize(trunc=trunc))
    end
    p
end

function build_K(c,lat,pos,trunc;den=false)
    if c.r==2
        return c.scale*polynomial_K(c,lat,pos,den ? c.k0masks : c.kmasks,den ? c.k0coeff : c.kcoeff,trunc)
    end
    p=coherent_K(c,lat,pos,den ? c.links0 : c.links,trunc)
    if c.method===:weyl
        for (a,bar) in pos[c.r+1:end];p=contact_pair(p,a,bar,trunc);end
    end
    c.scale*p
end

function coefficient_error(p,c,pos)
    isempty(c.kmasks)&&return NaN
    ms=[[Main.KernelAudit.site_matrix(p[j],s) for s in 0:1] for j in 1:length(p)]
    function coeff(m)
        inds=Int[pos[j][k+1] for j in eachindex(pos) for k in 0:1 if !iszero((m>>(2*(j-1)+k))&1)]
        s=(-1)^sum((inds[i]>inds[j] for i in eachindex(inds) for j in i+1:length(inds));init=0)
        v=ones(ComplexF64,1,1)
        for j in eachindex(ms);v=v*ms[j][j in inds ? 2 : 1];end
        s*v[1,1]*GTEMPO.scaling(p)^length(p)
    end
    maxerr=0.
    for (m,z) in zip(c.kmasks,c.kcoeff);maxerr=max(maxerr,abs(coeff(m)-c.scale*z));end
    maxerr
end

function run_one(c;cutoff=1e-13)
    tr=truncdimcutoff(D=160,ϵ=cutoff);lat,pos=layout(c)
    ti=@timed build_if(c,lat,pos,tr)
    tk=@timed build_K(c,lat,pos,tr)
    tz=@timed build_K(c,lat,pos,tr;den=true)
    tn=@timed integrate(lat,tk.value,ti.value)
    td=@timed integrate(lat,tz.value,ti.value)
    num=c.scalar*tn.value;den=c.scalar*td.value
    err=max(abs(num-c.num),abs(den-c.den))/max(1,abs(c.num),abs(c.den))
    ke=coefficient_error(tk.value,c,pos)
    println(c.name," same-target error=",err," Kcoeff=",ke);flush(stdout)
    @assert err<1e-9 (c.name,num,c.num,den,c.den)
    @assert isnan(ke)||ke<1e-10 (c.name,ke)
    (;name=c.name,N=c.N,r=c.r,axis=c.axis,method=c.method,green=c.green,
      same_target_error=err,K_coefficient_error=ke,physical_error=abs(num/den-c.physical),
      if_seconds=ti.time,Knum_seconds=tk.time,Kden_seconds=tz.time,contraction_seconds=tn.time+td.time,
      IFbond=bonddim(ti.value),Kbond=bonddim(tk.value),num_re=real(num),num_im=imag(num),den_re=real(den),den_im=imag(den))
end

function run()
    rows=[]
    open(joinpath(ROOT,"results","package.csv"),"w") do io
        for c in CASES
            r=run_one(c);isempty(rows)&&println(io,join(keys(r),','));push!(rows,r)
            println(io,join(values(r),','));flush(io)
        end
    end
    rows
end
end
