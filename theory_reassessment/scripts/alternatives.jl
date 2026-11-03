module TheoryAlternatives
using GTEMPO,QuAPI,LinearAlgebra
const ROOT=dirname(@__DIR__)
include(joinpath(ROOT,"results","alternative_cases.jl"))
function run()
    open(joinpath(ROOT,"results","alternatives.csv"),"w") do io
        println(io,"axis,N,method,g_error,Z_error,physical_error,bond")
        for c in CASES
            N=c.N;T=c.axis===:real ? .4 : .2
            lat=c.axis===:real ? GrassmannLattice(N=N,δt=T/N,contour=:real) : GrassmannLattice(Nt=N,Nτ=N,δt=T/N,δτ=1/N,contour=:mixed)
            tr=truncdimcutoff(D=160,ϵ=1e-13)
            p=vacuumstate(ComplexF64,lat)
            for i in axes(c.D,1)
                ids=[j for j in axes(c.D,2) if abs(c.D[i,j])>1e-16]
                isempty(ids)&&continue
                b,k=c.nodes[i];row=lat[GrassmannContourIndex(k,b,1,true)]
                cols=[lat[GrassmannContourIndex(c.nodes[j][2],c.nodes[j][1],1,false)] for j in ids]
                f=GTEMPO.partialmpo(row,cols,c.D[i,ids])*vacuumstate(ComplexF64,lat)
                p=mult!(p,f;trunc=tr)[1]
            end
            model=ToulouseIM(ϵ_d=.7)
            K=boundarycondition!(sysdynamics(lat,model;trunc=tr),lat)
            c.axis===:real&&(K=systhermalstate!(K,lat,model;trunc=tr,β=1.))
            cache=environments(lat,K,p)
            g=[cached_greater(lat,k,K,p;cache) for k in 1:lat.kt]
            Z=integrate(lat,K,p);eg=maximum(abs.(g-c.g));ez=abs(Z-c.Z)
            @assert max(eg,ez)<1e-9 (c.axis,N,c.method,eg,ez)
            ref=im*Main.ShiftStudy.reference(c.axis,N)
            println(io,join((c.axis,N,c.method,eg,ez,norm(g-ref)/norm(ref),bonddim(p)),','));flush(io)
        end
    end
    println("Native alternative contractions passed")
end
end
