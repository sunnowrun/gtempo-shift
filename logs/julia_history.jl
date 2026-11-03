# time: 2026-11-03 20:43:19 CST
# mode: julia
	using Pkg; Pkg.develop([PackageSpec(path=joinpath(pwd(), "vendor", p)) for p in ["ImpurityModelBase", "QuAPI", "ExpExp", "FiniteMPSAlgorithms", "Z2Tensors", "GTEMPO"]]); Pkg.instantiate(); println("ENV_READY ", getpid(), " ", gethostname())
# time: 2026-11-03 20:47:32 CST
# mode: julia
	Pkg.offline(true); Pkg.develop([PackageSpec(path=joinpath(pwd(),"vendor",p)) for p in ["ImpurityModelBase","QuAPI","ExpExp","FiniteMPSAlgorithms","Z2Tensors","GTEMPO"]]); Pkg.instantiate(); println("ENV_READY")
# time: 2026-11-03 20:48:05 CST
# mode: julia
	open("logs/study.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; include("scripts/study.jl"); r=ShiftStudy.case(:real,4,:current); println("SMOKE ",ShiftStudy.rel(r.g,r.reference)," ",r.g); ShiftStudy.sequence_test(); end; end; println("SMOKE_DONE")
# time: 2026-11-03 20:55:37 CST
# mode: julia
	open("logs/smoke.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; r=ShiftStudy.case(:real,4,:current); println("SMOKE ",ShiftStudy.rel(r.g,r.reference)," ",r.g); ShiftStudy.sequence_test(); end; end; println("SMOKE_DONE")
# time: 2026-11-03 20:59:02 CST
# mode: julia
	include("scripts/kernel_audit.jl"); include("scripts/link_if.jl"); include("scripts/delay_mpo.jl"); println("AUDIT_MODULES_READY")
# time: 2026-11-03 20:59:12 CST
# mode: julia
	open("logs/kernel.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; KernelAudit.coupling_scaling(); DelayMPO.run(); end; end; println("KERNEL_DONE"); include("scripts/run_convergence.jl"); println("CONVERGENCE_DONE")
# time: 2026-11-03 21:04:05 CST
# mode: julia
	include("scripts/repaired_tti.jl"); println("REPAIR_MODULE_READY")
# time: 2026-11-03 21:04:16 CST
# mode: julia
	open("logs/tti.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; KernelAudit.tti_audit(); end; end; println("TTI_AUDIT_DONE"); open("logs/link.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; LinkIF.run(); end; end; println("LINK_DONE")
# time: 2026-11-03 21:08:22 CST
# mode: julia
	include("scripts/link_if.jl"); RepairedTTI.install(); println("PATCH_READY")
# time: 2026-11-03 21:09:44 CST
# mode: julia
	open("logs/repaired.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; RepairedTTI.run(); end; end; RepairedTTI.restore(); println("REPAIRED_DONE"); open("logs/link_matched.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; LinkIF.run(); end; end; println("LINK_DONE")
# time: 2026-11-03 21:13:02 CST
# mode: julia
	include("scripts/repaired_tti.jl"); include("scripts/link_if.jl"); RepairedTTI.install(); println("PATCH2_READY")
# time: 2026-11-03 21:13:28 CST
# mode: julia
	include("scripts/repaired_tti.jl"); RepairedTTI.install(); println("PATCH2_READY")
# time: 2026-11-03 21:13:38 CST
# mode: julia
	open("logs/repaired2.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; RepairedTTI.run(); end; end; RepairedTTI.restore(); println("REPAIRED2_DONE"); open("logs/strang.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; LinkIF.run_strang(); end; end; println("STRANG_DONE")
# time: 2026-11-03 21:15:55 CST
# mode: julia
	open("logs/repaired2_error.log","w") do io; try; redirect_stdio(stdout=io,stderr=io) do; RepairedTTI.run(); end; catch err; showerror(io,err,catch_backtrace()); end; end; RepairedTTI.restore(); println("RESTORED"); open("logs/strang.log","w") do io; redirect_stdio(stdout=io,stderr=io) do; LinkIF.run_strang(); end; end; println("STRANG_DONE")
# time: 2026-11-03 21:17:06 CST
# mode: julia
	include("scripts/repaired_tti.jl"); RepairedTTI.install(); println("PATCH3_READY")
# time: 2026-11-03 21:17:19 CST
# mode: julia
	open("logs/repaired3.log","w") do io; try; redirect_stdio(stdout=io,stderr=io) do; RepairedTTI.run(); end; catch err; showerror(io,err,catch_backtrace()); finally; RepairedTTI.restore(); end; end; println("REPAIRED3_DONE")
# time: 2026-11-03 21:19:28 CST
# mode: julia
	include("scripts/repaired_tti.jl"); RepairedTTI.install(); println("FINAL_PATCH_READY")
# time: 2026-11-03 21:19:37 CST
# mode: julia
	open("logs/repaired_final.log","w") do io; try; redirect_stdio(stdout=io,stderr=io) do; RepairedTTI.run(); end; catch err; showerror(io,err,catch_backtrace()); finally; RepairedTTI.restore(); end; end; println("FINAL_REPAIR_DONE")
# time: 2026-11-03 21:20:57 CST
# mode: julia
	include("scripts/repaired_tti.jl"); RepairedTTI.install(); println("VALIDATION_READY")
# time: 2026-11-03 21:21:09 CST
# mode: julia
	open("logs/repaired_verified.log","w") do io; try; redirect_stdio(stdout=io,stderr=io) do; RepairedTTI.run(); end; catch err; showerror(io,err,catch_backtrace()); finally; RepairedTTI.restore(); end; end; println("VERIFIED_REPAIR_DONE")
