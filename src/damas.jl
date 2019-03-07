import Base: getindex, iterate

"""
    damas!(x, env::Environment, b [,f<:AbstractArray, ω::Real]; maxiter=10) -> x
Performs exactly `maxiter` DAMAS iterations for all frequency bins in `env.fn` or frequencies in `f`. Successive Over Relaxation (SOR) can be set with relaxation parameter `ω`.
Default is `ω=1` corresponding to no relaxation.
"""
function damas!(x, env, b, f::T=env.fn, ω::Real=1.0; maxiter::Int=10) where T <: AbstractArray
    # TODO: Check input sizes and index correct from f.
    for i = 1:env.Nf
        xf = x[:,i]
        iterable = DAMASSORIterable(xf, env.steeringvec.arr[:,:,i], zeros(env.N), zeros(env.N), b[:,i], ω, maxiter)
        for _ = iterable end
        x[:,i] = xf
    end
    x
end

function damas!(x, env, b, f::Number, ω::Real=1.0; maxiter::Int=10)
    fin = argmin(abs.(env.fn .- f))
    steer = env.steeringvec.arr[:,:,fin]
    println("Computing DAMAS for f = $(env.fn[fin])")
    iterable = DAMASSORIterable(x, steer, similar(x), similar(x), b[:,fin], ω, maxiter)
    for _ = iterable end
    x
end

mutable struct DAMASSORIterable{solT,steerT,vecT,rhsT,numT}
    x::solT
    steer::steerT
    tmp1::vecT
    tmp2::vecT
    b::rhsT
    ω::numT
    maxiter::Int
end

start(::DAMASSORIterable) = 1
done(it::DAMASSORIterable, iteration::Int) = iteration > it.maxiter
function iterate(s::DAMASSORIterable, iteration::Int=start(s))
    if done(s, iteration) return nothing end

    n = size(s.b, 1)

    for col = 1 : n
        AeroAcoustics.psf_col!(s.tmp2,s.steer,col)
        @simd for row = 1 : col - 1
            @inbounds s.tmp1[row] = max(0.0,s.tmp1[row] - s.tmp2[row] * s.x[col])
        end

        s.tmp1[col] = s.b[col]
    end

    for col = 1 : n
        AeroAcoustics.psf_col!(s.tmp2,s.steer,col)
        @inbounds s.x[col] += s.ω * (s.tmp1[col] / s.tmp2[col] - s.x[col])
        @simd for row = col + 1 : n
            @inbounds s.tmp1[row] = max(0.0,s.tmp1[row] - s.tmp2[row] * s.x[col])
        end
    end

    nothing, iteration + 1
end
