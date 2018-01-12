function fistaprox!(x::Array{T,2}, psf::Array{T,2}, b::Array{T,2}, g; tol::R=1e-3, maxit::Int=5000,lam::R=0.01) where {T <: AbstractFloat, R <: Real}
    x_prev = copy(x)
    const Nx,Ny = size(x)
    center = round(Int64,Nx/2)+1,round(Int64,Ny/2)+1
    Fps = fft(circshift(psf,center))
    FpsT = fft(circshift(psf',center))
    s = rand(size(psf))
    for k = 1:10
        s = ifft(fft(s).*Fps)/vecnorm(s)
    end
    const L = vecnorm(s)^2
    const beta = 1./L
    lambda = lam*norm(real(ifft(FpsT.*fft(b))),Inf)
    f = g(lambda)
    for it = 1:maxit
        x_extr = x + (it-2)/(it+1)*(x - x_prev)          # extrapolation step
        res = real(ifft(Fps.*fft(x_extr)))-b             # res = A*x_extr - b
        y = real(x_extr-beta*real(ifft(FpsT.*fft(res)))) # y = x_extr - gam*(A'*res)
        x_prev .= x                                      # store current iterate
        prox!(x, f, y, beta)
        nonneg!(x)
        #prox!(x, IndNonnegative(), y, beta)
        # stopping criterion
        if norm(x_extr-x, Inf)/beta <= tol*(1+norm(x, Inf))
            break
        end
    end
    return x
end

# For multidimensional arrays, this is per definition warm-starting:
function fistaprox!(x::Array{T,3}, psf::Array{T,3}, b::Array{T,3}, g; tol::R=1e-3, maxit::Int=5000,lam::R=0.01) where {T <: AbstractFloat, R <: Real}
    const Nx,Ny,Nf = size(b)
    y = zeros(Nx,Ny)
    for i in Nf:-1:1
        fistaprox!(y,psf[:,:,i],b[:,:,i],g;tol=tol, maxit=maxit,lam=lam)
        x[:,:,i] = y
    end
    return x
end

function nonneg!(x::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    for k in eachindex(x)
        if x[k] < 0
            x[k] = zero(T)
        else
            x[k] = x[k]
        end
    end
    return 0.0
end

# squared L2 norm (times a constant, or weighted)
#=
"""
**Squared Euclidean norm (weighted)**
    SqrNormL2(λ=1.0)
With a nonnegative scalar `λ`, returns the function
```math
f(x) = \\tfrac{λ}{2}\\|x\\|^2.
```
With a nonnegative array `λ`, returns the function
```math
f(x) = \\tfrac{1}{2}∑_i λ_i x_i^2.
```
"""

mutable struct SqrNormL2{T <: Union{Real, AbstractArray}} <: ProximableFunction
  lambda::T
  function SqrNormL2{T}(lambda::T) where {T <: Union{Real,AbstractArray}}
    if any(lambda .< 0)
      error("coefficients in λ must be nonnegative")
    else
      new(lambda)
    end
  end
end

is_convex(f::SqrNormL2) = true
is_smooth(f::SqrNormL2) = true
is_separable(f::SqrNormL2) = true
is_quadratic(f::SqrNormL2) = true
is_strongly_convex(f::SqrNormL2) = all(f.lambda .> 0)

SqrNormL2{T <: Real}(lambda::T=1.0) = SqrNormL2{T}(lambda)

SqrNormL2{T <: AbstractArray}(lambda::T) = SqrNormL2{T}(lambda)

function (f::SqrNormL2{S}){S <: Real, T <: Real}(x::AbstractArray{T})
  return (f.lambda/2)*vecnorm(x)^2
end

function (f::SqrNormL2{S}){S <: AbstractArray, T <: Real}(x::AbstractArray{T})
  sqnorm = 0.0
  for k in eachindex(x)
    sqnorm += f.lambda[k]*abs2(x[k])
  end
  return 0.5*sqnorm
end

prox!{S <: Real, T <: Real}(y::AbstractArray{T}, f::AeroAcoustics.SqrNormL2{S}, x::AbstractArray{T}, gamma::Real=1.0) = prox!(y, f::ProximalOperators.SqrNormL2{S}, x, gamma)

prox!{S <: AbstractArray, T <: Real}(y::AbstractArray{T}, f::AeroAcoustics.SqrNormL2{S}, x::AbstractArray{T}, gamma=1.0) = prox!(y, f::ProximalOperators.SqrNormL2{S}, x, gamma)

prox!{S <: Real, T <: Real}(y::AbstractArray{T}, f::AeroAcoustics.SqrNormL2{S}, x::AbstractArray{T}, gamma::AbstractArray) = prox!(y, f::ProximalOperators.SqrNormL2{S}, x, gamma)

prox!{S <: AbstractArray, T <: Real}(y::AbstractArray{T}, f::AeroAcoustics.SqrNormL2{S}, x::AbstractArray{T}, gamma::AbstractArray) = prox!(y, f::ProximalOperators.SqrNormL2{S}, x, gamma)

fun_name(f::SqrNormL2) = "weighted squared Euclidean norm"
fun_dom(f::SqrNormL2) = "AbstractArray{Real}, AbstractArray{Complex}"
fun_expr{T <: Real}(f::SqrNormL2{T}) = "x ↦ (λ/2)||x||^2"
fun_expr{T <: AbstractArray}(f::SqrNormL2{T}) = "x ↦ (1/2)sum( λ_i (x_i)^2 )"
fun_params{T <: Real}(f::SqrNormL2{T}) = "λ = $(f.lambda)"
fun_params{T <: AbstractArray}(f::SqrNormL2{T}) = string("λ = ", typeof(f.lambda), " of size ", size(f.lambda))

"""
**``L_1`` norm**
    NormL1(λ=1.0)
With a nonnegative scalar parameter λ, returns the function
```math
f(x) = λ\\cdot∑_i|x_i|.
```
With a nonnegative array parameter λ, returns the function
```math
f(x) = ∑_i λ_i|x_i|.
```
"""

mutable struct NormL1{T <: Union{Real, AbstractArray}} <: ProximableFunction
  lambda::T
  function NormL1{T}(lambda::T) where {T <: Union{Real, AbstractArray}}
    if !(eltype(lambda) <: Real)
      error("λ must be real")
    end
    if any(lambda .< 0)
      error("λ must be nonnegative")
    else
      new(lambda)
    end
  end
end

NormL1{R <: Real}(lambda::R=1.0) = NormL1{R}(lambda)

NormL1{A <: AbstractArray}(lambda::A) = NormL1{A}(lambda)

function (f::NormL1{R}){R <: Real}(x::AbstractArray)
  return f.lambda*vecnorm(x,1)
end

function (f::NormL1{A}){A <: AbstractArray}(x::AbstractArray)
  return vecnorm(f.lambda.*x,1)
end

prox!{A <: AbstractArray, R <: Real}(y::AbstractArray{R}, f::NormL1{A}, x::AbstractArray{R}, gamma::Real=1.0)

prox!{A <: AbstractArray, R <: Real}(y::AbstractArray{Complex{R}}, f::NormL1{A}, x::AbstractArray{Complex{R}}, gamma::Real=1.0)

prox!{T <: Real, R <: Real}(y::AbstractArray{R}, f::NormL1{T}, x::AbstractArray{R}, gamma::Real=1.0)

prox!{T <: Real, R <: Real}(y::AbstractArray{Complex{R}}, f::NormL1{T}, x::AbstractArray{Complex{R}}, gamma::Real=1.0)

prox!{A <: AbstractArray, R <: Real}(y::AbstractArray{R}, f::NormL1{A}, x::AbstractArray{R}, gamma::AbstractArray)

prox!{A <: AbstractArray, R <: Real}(y::AbstractArray{Complex{R}}, f::NormL1{A}, x::AbstractArray{Complex{R}}, gamma::AbstractArray)

prox!{T <: Real, R <: Real}(y::AbstractArray{R}, f::NormL1{T}, x::AbstractArray{R}, gamma::AbstractArray)

prox!{T <: Real, R <: Real}(y::AbstractArray{Complex{R}}, f::NormL1{T}, x::AbstractArray{Complex{R}}, gamma::AbstractArray)

fun_name(f::NormL1) = "weighted L1 norm"
fun_dom(f::NormL1) = "AbstractArray{Real}, AbstractArray{Complex}"
fun_expr{R <: Real}(f::NormL1{R}) = "x ↦ λ||x||_1"
fun_expr{A <: AbstractArray}(f::NormL1{A}) = "x ↦ sum( λ_i |x_i| )"
fun_params{R <: Real}(f::NormL1{R}) = "λ = $(f.lambda)"
fun_params{A <: AbstractArray}(f::NormL1{A}) = string("λ = ", typeof(f.lambda), " of size ", size(f.lambda))

# TODO Maybe not necessary...
# Extent ProximalOperators to include L1_pos

export NormL1Pos

mutable struct NormL1Pos{T <: Union{Real, AbstractArray}} <: ProximableFunction
    lambda::T
    function NormL1Pos{T}(lambda::T) where {T <: Union{Real, AbstractArray}}
        if !(eltype(lambda) <: Real)
            error("λ must be real")
        end
        if any(lambda .< 0)
            error("λ must be nonnegative")
        else
            new(lambda)
        end
    end
end

is_separable(f::NormL1Pos) = true
is_convex(f::NormL1Pos) = true


NormL1Pos{R <: Real}(lambda::R=1.0) = NormL1Pos{R}(lambda)
NormL1Pos{A <: AbstractArray}(lambda::A) = NormL1Pos{A}(lambda)

function (f::NormL1Pos{R}){R <: Real}(x::AbstractArray)
    return f.lambda*vecnorm(x,1)
end

function (f::NormL1Pos{A}){A <: AbstractArray}(x::AbstractArray)
    return vecnorm(f.lambda.*x,1)
end

function prox!(y::AbstractArray{R,N}, f::NormL1Pos{T}, x::AbstractArray{R,N}, gamma::Real=1.0) where {T <: Real, R <: Real, N}
    for k in eachindex(x)
        if x[k]-f.lambda*gamma < 0
            y[k] = zero(R)
        else
            y[k] = x[k]-f.lambda*gamma
        end
    end
    return 0.0
end

prox!(y::AbstractArray{R}, f::NormL1Pos, x::AbstractArray{R}, gamma::AbstractArray{R}) where R <: Real = prox!(y, f, x, 1.0)

fun_name(f::NormL1Pos) = "Non-negative weighted L1 norm"
fun_params{R <: Real}(f::NormL1Pos{R}) = "λ = $(f.lambda)"
fun_params{A <: AbstractArray}(f::NormL1Pos{A}) = string("λ = ", typeof(f.lambda), " of size ", size(f.lambda))
=#
