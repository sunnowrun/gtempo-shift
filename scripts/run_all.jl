# Start one Julia REPL from the workspace with --project=env, then:
# include("scripts/run_all.jl")
# All modules are loaded at top level so world ages advance between steps.
include("study.jl")
include("kernel_audit.jl")
include("link_if.jl")
include("delay_mpo.jl")
include("repaired_tti.jl")

function logged(f, name)
    open(joinpath(ShiftStudy.ROOT,"logs",name),"w") do io
        redirect_stdio(f; stdout=io,stderr=io)
    end
end

logged("sequence.log") do
    ShiftStudy.sequence_test()
    DelayMPO.run()
    KernelAudit.coupling_scaling()
end
logged(ShiftStudy.convergence,"convergence.log")
logged(KernelAudit.tti_audit,"tti.log")
logged(LinkIF.run,"link_matched.log")
logged(LinkIF.run_strang,"strang.log")

RepairedTTI.install()
try
    logged(RepairedTTI.run,"repaired2.log")
finally
    RepairedTTI.restore()
    ShiftStudy.MODE[]=:current
end
println("Study finished; run scripts/summarize.py on the CPU allocation to export the figure.")
