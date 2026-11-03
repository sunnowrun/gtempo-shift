module RepairedTTI
using GTEMPO, QuAPI, ImpurityModelBase, Z2Tensors, ExpExp, LinearAlgebra, Random
const ROOT=dirname(@__DIR__)
const SOURCE=read(joinpath(ROOT,"vendor","GTEMPO","src","influencefunctional","ttiif","util.jl"),String)
const START=findfirst("function _vertex_shifted_parts",SOURCE).start
const STOP=findnext("\nfunction _ti_mpotensor_shifted(",SOURCE,START).start
const ORIGINAL=SOURCE[START:prevind(SOURCE,STOP)]

function install()
    # Research patch, in memory only. The scalar sequence is shifted exactly;
    # samples on legs that will be dropped are excluded from the fit.
    # Preserve the raw generator's endpoint orientations. Its m1 family
    # generates upper lags, and its m2 family lower lags; the upstream
    # shifted implementation feeds those two families the opposite tails.
    cut=findfirst("\n\tif kind === :τ",ORIGINAL).start
    s=ORIGINAL[1:prevind(ORIGINAL,cut)] * "\n" * """
        if kind === :τ || kind === :+
            m1s=expand_decayterm(GenericDecayTerm(a1,abar1,-corr.ηₖⱼ[3:M],middle=JW4),alg=alg)
            m2s=expand_decayterm(GenericDecayTerm(abar2,a2,corr.ηⱼₖ[2:M-1],middle=JW4),alg=alg)
            m2s=[ExponentialDecayTerm(t.a,t.b;middle=t.m,α=t.α,coeff=t.coeff/t.α) for t in m2s]
            return corr.ηₖⱼ[2],m1s,m2s
        else
            m1s=expand_decayterm(GenericDecayTerm(a1,abar1,-corr.ηₖⱼ[2:M-1],middle=JW4),alg=alg)
            m1s=[ExponentialDecayTerm(t.a,t.b;middle=t.m,α=t.α,coeff=t.coeff/t.α) for t in m1s]
            m2s=expand_decayterm(GenericDecayTerm(abar2,a2,corr.ηⱼₖ[3:M],middle=JW4),alg=alg)
            return corr.ηⱼₖ[2],m1s,m2s
        end
    end
    """
    Core.eval(GTEMPO,Meta.parse(s))
end
restore()=Core.eval(GTEMPO,Meta.parse(ORIGINAL))

function sampled_coeff_error(psi,ref,lat,bs)
    function matrices(p)
        [ [Main.KernelAudit.site_matrix(p[j],σ) for σ in 0:1] for j in 1:length(p)]
    end
    function value(ms,occ)
        v=ones(ComplexF64,1,1)
        for j in eachindex(ms);v=v*ms[j][j in occ ? 2 : 1];end
        v[1,1]
    end
    a,b=matrices(psi),matrices(ref)
    ca,cb=value(a,Int[]),value(b,Int[])
    bars=[lat[l] for l in bs if l.conj];anns=[lat[l] for l in bs if !l.conj]
    rng=MersenneTwister(90210)
    maxerr=0.0
    for q in 1:256
        rank=1+mod(q,4)
        occ=vcat(shuffle(rng,bars)[1:rank],shuffle(rng,anns)[1:rank])
        maxerr=max(maxerr,abs(value(a,occ)/ca-value(b,occ)/cb))
    end
    maxerr
end

function run()
    Main.ShiftStudy.MODE[]=:current
    tr=truncdimcutoff(D=160,ϵ=1e-13)
    alg=ExactTTIIF(alg_expan=OverDeterminedProny(n=20,tol=1e-12),alg_mult=SVDCompression(trunc=tr))
    open(joinpath(ROOT,"results","repaired_tti.csv"),"w") do io
        println(io,"N,stage,relative_kernel_error,max_abs_error,relative_if_error,off_contact_max_error,sampled_if_coeff_max_error")
        for N in (3,4,8)
            lat=GrassmannLattice(N=N,δt=0.4/N,contour=:real)
            corr=correlationfunction(Main.ShiftStudy.bath(),lat)
            bs,C=Main.KernelAudit.accessor_matrix(corr)
            ref=hybriddynamics(lat,corr,trunc=tr)
            psi=Base.invokelatest(hybriddynamics,lat,corr,alg)
            D=Main.KernelAudit.quadratic_matrix(psi,lat,bs)
            contacts=Tuple{Int,Int}[]
            for (i,r) in enumerate(bs),(j,c) in enumerate(bs)
                r.conj && !c.conj || continue
                localterm = r.branch===c.branch &&
                    (r.j==c.j || (r.branch===:+ && r.j==c.j+1) || (r.branch===:- && c.j==r.j+1))
                window = (r.branch===:+ && r.j==1) || (r.branch===:- && r.j==lat.k) ||
                    (c.branch===:+ && c.j==lat.k) || (c.branch===:- && c.j==1)
                (localterm || window) && push!(contacts,(i,j))
            end
            off=copy(D-C)
            for (i,j) in contacts; off[i,j]=0; end
            println(io,join((N,"shifted_tail",norm(D-C)/norm(C),maximum(abs.(D-C)),norm(psi-ref)/norm(ref),maximum(abs.(off)),sampled_coeff_error(psi,ref,lat,bs)),','));flush(io)
            # Diagnostic extraction of local and boundary coefficients.
            # Production coefficients can be read directly from the fitted
            # poles; extracting the MPS kernel here is an independent check.
            for (i,j) in contacts
                r,c=bs[i],bs[j]
                pr=lat[GrassmannContourIndex(r.j,r.branch,1,true)]
                pc=lat[GrassmannContourIndex(c.j,c.branch,1,false)]
                apply!(exp(GTerm(pr,pc;coeff=C[i,j]-D[i,j])),psi)
                canonicalize!(psi,alg=GTEMPO.Orthogonalize(trunc=tr))
            end
            D=Main.KernelAudit.quadratic_matrix(psi,lat,bs)
            e=norm(D-C)/norm(C)
            ei=norm(psi-ref)/norm(ref)
            es=sampled_coeff_error(psi,ref,lat,bs)
            println(io,join((N,"tail_plus_exact_local_and_boundary",e,maximum(abs.(D-C)),ei,maximum(abs.(off)),es),','));flush(io)
            @assert e<1e-8
            @assert es<1e-10
            println("repaired TTI N=",N," kernel error=",e," IF error=",ei);flush(stdout)
        end
    end
end
end
