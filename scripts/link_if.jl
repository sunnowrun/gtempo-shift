module LinkIF
using GTEMPO, QuAPI, ImpurityModelBase, Z2Tensors, LinearAlgebra
const ROOT=dirname(@__DIR__)

# A reference constructor. It uses the existing public exact-cell/Schur API,
# then analytic row exponentials; no arbitrary kernel-index shifts are used.
# It is deliberately not advertised as an efficient infinite-bath TTI solver.
function build_mps(lat, hmat; trunc=Main.ShiftStudy.TRUNC)
    psi=vacuumstate(ComplexF64,lat)
    bs=hmat.basis
    for (i,row) in enumerate(bs)
        row.conj || continue
        ids=[j for j in eachindex(bs) if !bs[j].conj && !iszero(hmat.Δ[i,j])]
        isempty(ids) && continue
        index(l)=lat[GrassmannContourIndex(l.branch===:τ ? l.j+1 : l.j,l.branch,l.band,l.conj)]
        mpo=GTEMPO.partialmpo(index(row),[index(bs[j]) for j in ids],hmat.Δ[i,ids])
        factor=mpo*vacuumstate(ComplexF64,lat)
        psi=mult!(psi,factor,trunc=trunc)[1]
    end
    return psi
end

function real_case(N; ed=0.7, beta=1.0, T=0.4, weight=0.5, omega=1.0, matched=true, scheme=:lie)
    lat=GrassmannLattice(N=N,δt=T/N,contour=:real)
    spec=DiracDelta(ω=omega,α=weight)
    ex=exact_fermionic_Δt(spec;β=beta,t=T,N=N)
    C=copy(ex.hmat.Δ)
    if matched
        x=exp(-im*T/N*ed)
        for (i,leg) in enumerate(ex.hmat.basis)
            # B*A: A acts on the incoming impurity leg on +, and on
            # the outgoing conjugate leg on - for A†*B†.
            if scheme===:strang
                half=exp(-im*T/N*ed/2)
                if (leg.branch===:+ && leg.conj && leg.j>=2)
                    C[i,:].*=half
                elseif (leg.branch===:+ && !leg.conj && leg.j<=N)
                    C[:,i].*=half
                elseif (leg.branch===:- && leg.conj && leg.j<=N)
                    C[i,:].*=conj(half)
                elseif (leg.branch===:- && !leg.conj && leg.j>=2)
                    C[:,i].*=conj(half)
                end
            elseif leg.j<=N && leg.branch===:+ && !leg.conj
                C[:,i].*=x
            elseif leg.j<=N && leg.branch===:- && leg.conj
                C[i,:].*=conj(x)
            end
        end
    end
    hm=QuAPI.Dense.GrassmannHybridizationMatrix(ex.hmat.basis,C)
    psiI=build_mps(lat,hm)
    model=ToulouseIM(ϵ_d=ed)
    psiK=boundarycondition!(sysdynamics(lat,model,trunc=Main.ShiftStudy.TRUNC),lat,trunc=Main.ShiftStudy.TRUNC)
    psiK=systhermalstate!(psiK,lat,model,β=beta,trunc=Main.ShiftStudy.TRUNC)
    cache=environments(lat,psiK,psiI)
    g=[-im*cached_greater(lat,k,psiK,psiI,cache=cache) for k in 1:lat.k]
    h=[ed sqrt(weight);sqrt(weight) omega]
    n0=1/(1+exp(beta*ed))
    ref=[-im*exp(-im*t*h)[1,1]*(1-n0) for t in range(0,T,length=N+1)]
    # Microscopic split chain is an independent single-particle oracle.
    hb=[0.0 sqrt(weight);sqrt(weight) omega]
    hs=Diagonal([ed,0.0])
    U=scheme===:strang ? exp(-im*T/N/2*hs)*exp(-im*T/N*hb)*exp(-im*T/N/2*hs) :
        exp(-im*T/N*hb)*exp(-im*T/N*hs)
    split=[-im*(U^k)[1,1]*(1-n0) for k in 0:N]
    return (;g,ref,split,Z=integrate(lat,psiK,psiI)*ex.prefactor,bond=bonddim(psiI))
end

function run()
    open(joinpath(ROOT,"results","link_if.csv"),"w") do io
        println(io,"ed,N,relative_error,split_oracle_error,Z_re,Z_im,bond,seconds")
        for ed in (0.0,0.7), N in (4,8,16,32)
            println("exact-link case ed=",ed," N=",N);flush(stdout)
            sec=@elapsed r=real_case(N;ed=ed)
            e=norm(r.g-r.ref)/norm(r.ref)
            eo=norm(r.g-r.split)/norm(r.split)
            println(io,join((ed,N,e,eo,real(r.Z),imag(r.Z),r.bond,sec),','));flush(io)
            println("  error=",e," split mismatch=",eo," Z=",r.Z);flush(stdout)
            @assert eo<1e-8
        end
    end
    open(joinpath(ROOT,"results","link_if_naive_control.csv"),"w") do io
        println(io,"ed,N,relative_error,split_oracle_error,Z_re,Z_im,bond,seconds")
        for N in (4,8,16,32)
            sec=@elapsed r=real_case(N;ed=0.7,matched=false)
            println(io,join((0.7,N,norm(r.g-r.ref)/norm(r.ref),norm(r.g-r.split)/norm(r.split),real(r.Z),imag(r.Z),r.bond,sec),','));flush(io)
        end
    end
end

function run_strang()
    open(joinpath(ROOT,"results","link_if_strang.csv"),"w") do io
        println(io,"ed,N,relative_error,split_oracle_error,Z_re,Z_im,bond,seconds")
        for N in (4,8,16,32)
            sec=@elapsed r=real_case(N;ed=0.7,scheme=:strang)
            eo=norm(r.g-r.split)/norm(r.split)
            println(io,join((0.7,N,norm(r.g-r.ref)/norm(r.ref),eo,real(r.Z),imag(r.Z),r.bond,sec),','));flush(io)
            println("Strang N=",N," error=",norm(r.g-r.ref)/norm(r.ref)," split mismatch=",eo);flush(stdout)
            @assert eo<1e-8
        end
    end
end
end
