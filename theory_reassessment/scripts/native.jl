module TheoryNative
using GTEMPO,QuAPI,LinearAlgebra
const ROOT=dirname(@__DIR__)
function run()
    open(joinpath(ROOT,"results","native.csv"),"w") do io
        println(io,"axis,N,method,k,g_re,g_im,Z_re,Z_im,seconds,bond")
        for axis in (:real,:imag,:mixed),N in (2,4,8,16),mode in (:raw,:vertex)
            sec=@elapsed r=Main.ShiftStudy.case(axis,N,mode)
            for k in eachindex(r.g)
                # Main study g includes -i on real/mixed; Gaussian oracle is C.
                g=axis===:imag ? r.g[k] : im*r.g[k]
                println(io,join((axis,N,mode,k,real(g),imag(g),real(r.Z),imag(r.Z),sec,r.bond),','))
            end
            flush(io);println("native ",axis," ",N," ",mode," ",sec);flush(stdout)
        end
    end
    Main.ShiftStudy.MODE[]=:current
end
end
