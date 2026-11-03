module TheoryExport
using GTEMPO, QuAPI, ImpurityModelBase, LinearAlgebra, DelimitedFiles
const ROOT=dirname(@__DIR__)

# Explicit read of the upstream accessors, without overwriting package methods.
function coefficient(corr,b1,j1,b2,j2,axis,mode)
    if mode===:vertex
        m=axis===:real ? QuAPI._vs_rt_acc(b1,j1,b2,j2,size(branch(corr;b1=:+,b2=:+),1)-1) :
          axis===:imag ? QuAPI._vs_τ_acc(j1,j2,size(branch(corr;b1=:τ,b2=:τ),1)) :
          QuAPI._vs_mixed_acc(b1,j1,b2,j2,QuAPI.isize(corr),QuAPI.rsize(corr)-1)
        m===nothing && return 0.0im
        v=m[5]*branch(corr;b1=m[1],b2=m[3])[m[2],m[4]]
    else
        v=branch(corr;b1=b1,b2=b2)[j1,j2]
    end
    axis===:real && b1!==b2 ? -v : v
end

function export_case(axis,N;omega=1.0,T=axis===:mixed ? .2 : .4,beta=1.,ed=.7,weight=.5,Ni=N,tag="main")
    spec=DiracDelta(ω=omega,α=weight)
    corr=axis===:real ? fermionic_Δt(spec;β=beta,t=T,N=N) :
         axis===:imag ? fermionic_Δτ(spec;β=beta,N=Ni) :
         fermionic_Δm(spec;β=beta,t=T,Nt=N,Nτ=Ni)
    nodes=axis===:real ? [(b,j) for b in (:+,:-) for j in 1:N+1] :
          axis===:imag ? [(:τ,j) for j in 1:Ni+1] :
          vcat([(:τ,j) for j in 1:Ni+1],[(b,j) for b in (:+,:-) for j in 1:N+1])
    id="$(tag)_$(axis)_$(N)_$(Ni)"
    open(joinpath(ROOT,"results",id*"_nodes.csv"),"w") do io
        for (b,j) in nodes;println(io,"$(b),$(j)");end
    end
    for mode in (:raw,:vertex)
        D=zeros(ComplexF64,length(nodes),length(nodes))
        for (r,(b,j)) in enumerate(nodes),(c,(b2,k)) in enumerate(nodes)
            (b===:τ && j==1 || b2===:τ && k==1) && continue
            D[r,c]=coefficient(corr,b,b===:τ ? j-1 : j,b2,b2===:τ ? k-1 : k,axis,mode)
        end
        writedlm(joinpath(ROOT,"results",id*"_$(mode).csv"),hcat(real.(D),imag.(D)),',')
    end
    (;id,axis,N,Ni,T,beta,ed,omega,weight)
end

function run()
    rows=[]
    for axis in (:real,:imag,:mixed),N in (1,2,4,8,16,32,64,128)
        axis!==:real && N==1 && continue # old τ shift reads [1,2], undefined for Nτ=1
        push!(rows,export_case(axis,N))
    end
    for omega in (-1.,0.),axis in (:real,:imag,:mixed),N in (4,8,16,32,64)
        push!(rows,export_case(axis,N;omega,tag=omega==0 ? "zero" : "filled"))
    end
    for T in (.05,1e-6),N in (4,8,16,32,64)
        push!(rows,export_case(:mixed,N;T,tag=T<1e-5 ? "short0" : "short"))
    end
    open(joinpath(ROOT,"results","cases.csv"),"w") do io
        println(io,join(keys(first(rows)),','))
        for r in rows;println(io,join(values(r),','));end
    end
    println("exported ",length(rows)," upstream cases")
end
end
