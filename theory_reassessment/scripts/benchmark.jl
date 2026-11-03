module TheoryBenchmark
using Statistics
const ROOT=dirname(@__DIR__)
function run()
    open(joinpath(ROOT,"results","warm_timings.csv"),"w") do io
        println(io,"axis,N,method,median_total_seconds,median_IF_seconds,median_K_seconds,median_contract_seconds,IFbond,Kbond,repetitions")
        for axis in (:real,:mixed),method in (:q,:weyl)
            c=only(filter(c->c.axis===axis&&c.method===method&&c.r==1&&c.N==8&&!c.green,Main.TheoryPackage.CASES))
            Main.TheoryPackage.run_one(c) # warm this exact specialization
            rows=[Main.TheoryPackage.run_one(c) for _ in 1:3]
            totals=[r.if_seconds+r.Knum_seconds+r.Kden_seconds+r.contraction_seconds for r in rows]
            println(io,join((axis,8,method,median(totals),median(r.if_seconds for r in rows),median(r.Knum_seconds+r.Kden_seconds for r in rows),median(r.contraction_seconds for r in rows),maximum(r.IFbond for r in rows),maximum(r.Kbond for r in rows),3),','));flush(io)
        end
    end
end
end
