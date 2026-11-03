module KernelAudit
using GTEMPO, QuAPI, ImpurityModelBase, Z2Tensors, ExpExp, LinearAlgebra
const ROOT=dirname(@__DIR__)

function accessor_matrix(corr)
    bs=QuAPI.Dense.basis(corr)
    C=zeros(ComplexF64,length(bs),length(bs))
    rt=QuAPI.branches(corr)==(:+, :-)
    for (i,r) in enumerate(bs), (j,c) in enumerate(bs)
        r.conj && !c.conj || continue
        C[i,j]=corr[conj(r),conj(c)]*(rt && r.branch!==c.branch ? -1 : 1)
    end
    return bs,C
end

function coupling_scaling()
    open(joinpath(ROOT,"results","kernel_coupling_scaling.csv"),"w") do io
        println(io,"axis,method,V,absolute_kernel_error,exact_kernel_norm")
        for axis in (:imag,:real), mode in (:raw,:vertex), V in (0.4,0.2,0.1,0.05)
            Main.ShiftStudy.MODE[]=mode
            spec=DiracDelta(ω=1.0,α=V^2)
            corr=axis===:imag ? fermionic_Δτ(spec;β=1.0,N=8) : fermionic_Δt(spec;β=1.0,t=0.4,N=8)
            ex=axis===:imag ? exact_fermionic_Δτ(spec;β=1.0,N=8) : exact_fermionic_Δt(spec;β=1.0,t=0.4,N=8)
            bs,C=accessor_matrix(corr)
            @assert bs==ex.hmat.basis
            println(io,join((axis,mode,V,norm(C-ex.hmat.Δ),norm(ex.hmat.Δ)),',')); flush(io)
        end
    end
    Main.ShiftStudy.MODE[]=:current
end

# Read a coefficient directly from the MPS, without materializing 2^L entries.
function offset(V,c)
    o=0
    for d in Z2Tensors.sectors(V)
        d.n<c.n || break
        o+=Z2Tensors.dim(V,d)
    end
    return o
end
function site_matrix(A,σ)
    Vl,Vr=GTEMPO.space_l(A),GTEMPO.space_r(A)
    M=zeros(GTEMPO.scalartype(A),Z2Tensors.dim(Vl),Z2Tensors.dim(Vr))
    for (fl,fr) in Z2Tensors.fusiontrees(A)
        sl,sp=fl.uncoupled
        sp.n==σ || continue
        fl.coupled.n==fr.coupled.n || continue
        sr=fr.uncoupled[1]
        ol,orr=offset(Vl,sl),offset(Vr,sr)
        dl,dr=Z2Tensors.dim(Vl,sl),Z2Tensors.dim(Vr,sr)
        M[ol+1:ol+dl,orr+1:orr+dr].=reshape(collect(A[fl,fr]),dl,dr)
    end
    return M
end
function quadratic_matrix(psi,lat,bs)
    ms=[[site_matrix(psi[p],σ) for σ in 0:1] for p in 1:length(psi)]
    function coeff(occupied)
        v=ones(ComplexF64,1,1)
        for p in eachindex(ms)
            v=v*ms[p][p in occupied ? 2 : 1]
        end
        return v[1,1]
    end
    c0=coeff(Int[])
    C=zeros(ComplexF64,length(bs),length(bs))
    for (i,r) in enumerate(bs), (j,c) in enumerate(bs)
        r.conj && !c.conj || continue
        pr=lat[GrassmannContourIndex(r.j,r.branch,1,true)]
        pc=lat[GrassmannContourIndex(c.j,c.branch,1,false)]
        C[i,j]=(pr<pc ? 1 : -1)*coeff([pr,pc])/c0
    end
    return C
end

function tti_audit()
    Main.ShiftStudy.MODE[]=:current
    tr=truncdimcutoff(D=160,ϵ=1e-13)
    alg=ExactTTIIF(alg_expan=OverDeterminedProny(n=20,tol=1e-12),alg_mult=SVDCompression(trunc=tr))
    open(joinpath(ROOT,"results","tti_kernel_audit.csv"),"w") do io
        println(io,"N,method,relative_kernel_error,max_abs_error,if_bond")
        for N in (3,4,8)
            lat=GrassmannLattice(N=N,δt=0.4/N,contour=:real)
            corr=correlationfunction(Main.ShiftStudy.bath(),lat)
            bs,C=accessor_matrix(corr)
            partial=hybriddynamics(lat,corr,trunc=tr)
            tti=hybriddynamics(lat,corr,alg)
            for (name,psi) in (("PartialIF",partial),("ExactTTIIF",tti))
                D=quadratic_matrix(psi,lat,bs)
                println(io,join((N,name,norm(D-C)/norm(C),maximum(abs.(D-C)),bonddim(psi)),',')); flush(io)
            end
            if N==4
                open(joinpath(ROOT,"results","tti_matrix_entries.csv"),"w") do mio
                    println(mio,"row_branch,row_j,col_branch,col_j,target_re,target_im,tti_re,tti_im")
                    D=quadratic_matrix(tti,lat,bs)
                    for (i,r) in enumerate(bs), (j,c) in enumerate(bs)
                        r.conj && !c.conj || continue
                        println(mio,join((r.branch,r.j,c.branch,c.j,real(C[i,j]),imag(C[i,j]),real(D[i,j]),imag(D[i,j])),','))
                    end
                end
            end
        end
    end
end

end
