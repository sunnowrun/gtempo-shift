module ShiftStudy
using GTEMPO, QuAPI, ImpurityModelBase, Z2Tensors, ExpExp
using LinearAlgebra, Printf, Dates

const ROOT = dirname(@__DIR__)
const MODE = Ref(:current)
const TRUNC = truncdimcutoff(D=160, ϵ=1e-13)
const BETA = 1.0
const ED = 0.7
const OMEGA = 1.0
const WEIGHT = 0.5
const H = [ED sqrt(WEIGHT); sqrt(WEIGHT) OMEGA]
const FERMI = inv(I + exp(BETA * H))

# Only this research process changes the accessor. The vendor source stays intact.
@eval QuAPI function Base.getindex(x::AbstractCorrelationFunction,
        i::GrassmannContourIndex, j::GrassmannContourIndex)
    D = branch(x; b1=branch(i), b2=branch(j))
    i.conj == j.conj && return zero(eltype(D))
    mode = Main.ShiftStudy.MODE[]
    bs = branches(x)
    vertex = mode === :vertex || (mode === :current && bs == (:+, :-))
    if vertex
        m = bs == (:τ,) ? _vs_τ_acc(i.j, j.j, size(D, 1)) :
            _vs_rt_acc(branch(i), i.j, branch(j), j.j, size(D, 1)-1)
        m === nothing && return zero(eltype(D))
        v = m[5] * branch(x; b1=m[1], b2=m[3])[m[2], m[4]]
    else
        v = D[i.j, j.j]
    end
    return i.conj ? -v : v
end
@eval QuAPI function Base.getindex(x::AbstractMixedCorrelationFunction,
        i::GrassmannContourIndex, j::GrassmannContourIndex)
    D = branch(x; b1=branch(i), b2=branch(j))
    i.conj == j.conj && return zero(eltype(D))
    if Main.ShiftStudy.MODE[] === :vertex
        nτ = size(branch(x; b1=:τ, b2=:τ), 1)
        nt = size(branch(x; b1=:+, b2=:+), 1)-1
        m = _vs_mixed_acc(branch(i), i.j, branch(j), j.j, nτ, nt)
        m === nothing && return zero(eltype(D))
        v = m[5] * branch(x; b1=m[1], b2=m[3])[m[2], m[4]]
    else
        v = D[i.j, j.j]
    end
    return i.conj ? -v : v
end

bath() = fermionicbath(DiracDelta(ω=OMEGA, α=WEIGHT), β=BETA)
rel(a,b) = norm(a-b)/norm(b)

function reference(axis, N)
    T = axis === :mixed ? 0.2 : 0.4
    ts = collect(range(0, T, length=N+1))
    if axis === :imag
        return ComplexF64[(exp(-τ*H)*(I-FERMI))[1,1] for τ in range(0,BETA,length=N+1)]
    elseif axis === :mixed
        return ComplexF64[-im*(exp(-im*t*H)*(I-FERMI))[1,1] for t in ts]
    else
        n0 = 1/(1+exp(BETA*ED))
        return ComplexF64[-im*exp(-im*t*H)[1,1]*(1-n0) for t in ts]
    end
end

function case(axis, N, mode; alg=nothing, trunc=TRUNC)
    MODE[] = mode
    T = axis === :mixed ? 0.2 : 0.4
    lat = axis === :imag ? GrassmannLattice(N=N, δτ=BETA/N, contour=:imag) :
        axis === :real ? GrassmannLattice(N=N, δt=T/N, contour=:real) :
        GrassmannLattice(Nt=N, Nτ=N, δt=T/N, δτ=BETA/N, contour=:mixed)
    corr = correlationfunction(bath(), lat)
    psiI = isnothing(alg) ? hybriddynamics(lat,corr,trunc=trunc) : hybriddynamics(lat,corr,alg)
    model = ToulouseIM(ϵ_d=ED)
    psiK = boundarycondition!(sysdynamics(lat,model,trunc=trunc),lat)
    axis === :real && (psiK=systhermalstate!(psiK,lat,model,trunc=trunc,β=BETA))
    cache = environments(lat,psiK,psiI)
    if axis === :imag
        g = cached_gf_fast(lat,psiK,psiI;c1=false,c2=true,b1=:τ,b2=:τ,cache=cache)
    else
        g = [-im*cached_greater(lat,k,psiK,psiI,cache=cache) for k in 1:lat.kt]
    end
    return (; g, reference=reference(axis,N), Z=integrate(lat,psiK,psiI),
            bond=bonddim(psiI), lat, corr, psiI, psiK)
end

function convergence(; small=false)
    data = Dict()
    open(joinpath(ROOT,"results","convergence.csv"),"w") do io
        println(io,"axis,method,N,step,relative_error,Z_re,Z_im,bond,seconds")
        open(joinpath(ROOT,"results","greens.csv"),"w") do gio
            println(gio,"axis,method,N,k,g_re,g_im,ref_re,ref_im")
            for axis in (:real,:imag,:mixed), mode in (:raw,:vertex)
                ns = small ? [4,8] : axis === :imag ? [8,16,32,64] : [4,8,16,32]
                for N in ns
                    println(now()," case ",axis," ",mode," N=",N); flush(stdout)
                    seconds=@elapsed r=case(axis,N,mode)
                    step=(axis===:imag ? BETA : axis===:mixed ? 0.2 : 0.4)/N
                    e=rel(r.g,r.reference)
                    println(io,join((axis,mode,N,step,e,real(r.Z),imag(r.Z),r.bond,seconds),',')); flush(io)
                    for k in eachindex(r.g)
                        println(gio,join((axis,mode,N,k,real(r.g[k]),imag(r.g[k]),real(r.reference[k]),imag(r.reference[k])),','))
                    end
                    flush(gio)
                    data[(axis,mode,N)]=r.g
                    println("  relerr=",e," bond=",r.bond," seconds=",seconds); flush(stdout)
                    if haskey(data,(axis,mode,N÷2))
                        extrap=2 .* r.g[1:2:end] .- data[(axis,mode,N÷2)]
                        er=rel(extrap,reference(axis,N÷2))
                        println(io,join((axis,Symbol(mode,:_richardson),N÷2,2step,er,NaN,NaN,r.bond,seconds),',')); flush(io)
                    end
                end
            end
        end
    end
    MODE[]=:current
    return data
end

# Delayed exponential sequence, k(d)=E at d=1, k(d)=tail(d-1) at d>=2.
# Tail fits use the author's latest ExpExp; the contact never enters Prony.
function contact_fit(tail, E; tol=1e-12)
    c, λ = exponential_expansion(collect(tail),alg=OverDeterminedProny(n=20,tol=tol))
    firstweights = c .* λ
    return (; E, weights=firstweights, poles=λ,
        evaluate=d->d==1 ? E : sum(firstweights .* λ.^(d-2)),
        error=norm([sum(firstweights .* λ.^(d-1)) for d in eachindex(tail)]-tail)/norm(tail))
end

function sequence_test()
    open(joinpath(ROOT,"results","prony_sequences.csv"),"w") do io
        println(io,"spectrum,N,tail_fit_error,delayed_max_error,contact,naive_contact_error,channels")
        for spec in (DiracDelta(ω=1.0,α=0.5),
                     DiscreteSpectrum([-1.4,0.3,1.7],[0.1,0.2,0.15])), N in (8,16,32)
            corr=fermionic_Δt(spec;β=1.0,t=0.4,N=N)
            η=QuAPI.branch(corr;b1=:+,b2=:+)
            tail=η.ηⱼₖ[2:end-1]
            E=η[1,1]
            f=contact_fit(tail,E)
            target=vcat(E,tail)
            err=maximum(abs.([f.evaluate(d) for d in eachindex(target)]-target))
            naive=tail[1]
            println(io,join((spec isa DiracDelta ? "one_mode" : "three_modes",N,f.error,err,abs(E),abs(E-naive),length(f.poles)),','))
            @assert err < 1e-10
        end
    end
end

end # module
