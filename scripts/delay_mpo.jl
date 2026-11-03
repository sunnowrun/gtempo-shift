module DelayMPO
using LinearAlgebra, Random
const ROOT=dirname(@__DIR__)

# Finite-state MPO for E δ[d,1] + 1[d>=2] Σ w_l λ_l^(d-2).
# Scalar probe variables verify the full pair sum; in a Grassmann MPO,
# replace x,y by endpoint operators and intermediate 1 by the JW string.
function pair_sum(x,y,E,w,λ)
    r=length(λ)
    v=zeros(ComplexF64,r+3); v[1]=1
    for j in eachindex(x)
        W=zeros(ComplexF64,r+3,r+3)
        W[1,1]=W[end,end]=1
        W[1,2]=x[j]
        W[2,end]=E*y[j]
        for l in 1:r
            W[2,l+2]=w[l]
            W[l+2,l+2]=λ[l]
            W[l+2,end]=y[j]
        end
        v=transpose(W)*v
    end
    return v[end]
end

function run()
    rng=MersenneTwister(128)
    open(joinpath(ROOT,"results","delay_mpo.csv"),"w") do io
        println(io,"poles,N,absolute_pair_sum_error,max_kernel_error,bond")
        for λ in (ComplexF64[0.0],ComplexF64[0.8+0.1im],ComplexF64[0.0,0.6,0.9-0.15im]), N in (3,8,32,128)
            E=0.25-0.13im
            w=randn(rng,ComplexF64,length(λ))
            x=randn(rng,ComplexF64,N);y=randn(rng,ComplexF64,N)
            kernel(d)=d==1 ? E : sum(w.*λ.^(d-2))
            target=sum(x[i]*y[j]*kernel(j-i) for i in 1:N-1 for j in i+1:N)
            actual=pair_sum(x,y,E,w,λ)
            err=abs(target-actual)
            ek=0.0
            for d in 1:N-1
                xp=zeros(N);yp=zeros(N);xp[1]=yp[d+1]=1
                ek=max(ek,abs(pair_sum(xp,yp,E,w,λ)-kernel(d)))
            end
            println(io,join((length(λ),N,err,ek,length(λ)+3),','));flush(io)
            @assert err<1e-10
            @assert ek<1e-12
        end
    end
end
end
