using Pkg
Pkg.activate(joinpath(@__DIR__,"..","env"))
root=dirname(@__DIR__)
packages=["ImpurityModelBase","QuAPI","ExpExp","FiniteMPSAlgorithms","Z2Tensors","GTEMPO"]
Pkg.develop([PackageSpec(path=joinpath(root,"vendor",p)) for p in packages])
Pkg.instantiate()
