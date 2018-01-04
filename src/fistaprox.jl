"""
** Non-negative ``L_1`` norm**
```
x_fistaprox = zeros(Env.Nx,Env.Ny)
x_fistaproxL1 = zeros(Env.Nx,Env.Ny)
x_fistaproxL1Pos = zeros(Env.Nx,Env.Ny)
fistaprox!(x_fistaprox, PSF[:,:,10], b[:,:,10], IndNonnegative(); tol=1e-6, maxit=1000)
fistaprox!(x_fistaproxL1, PSF[:,:,10], b[:,:,10], NormL1(1e-2); tol=1e-6, maxit=1000)
fistaprox!(x_fistaproxL1, PSF[:,:,10], b[:,:,10], NormL1Pos(1e-2); tol=1e-6, maxit=1000)
```
"""

function fistaprox!(x::Array{T,2}, psf::Array{T,2}, b::Array{T,2}, g; tol=1e-3, maxit=5000,lam=0.01) where T <: AbstractFloat
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
    g.lambda = lam*norm(real(ifft(FpsT.*fft(b))),Inf)
    for it = 1:maxit
        x_extr = x + (it-2)/(it+1)*(x - x_prev)          # extrapolation step
        res = real(ifft(Fps.*fft(x_extr)))-b             # res = A*x_extr - b
        y = real(x_extr-beta*real(ifft(FpsT.*fft(res)))) # y = x_extr - gam*(A'*res)
        x_prev .= x                                      # store current iterate
        prox!(x, g, y, beta)
        # stopping criterion
        if norm(x_extr-x, Inf)/beta <= tol*(1+norm(x, Inf))
            break
        end
    end
    return x
end

function fistaprox!(x::Array{T,3}, psf::Array{T,3}, b::Array{T,3}, g; tol=1e-3, maxit=5000,lam=0.01) where T <: AbstractFloat
    const Nx,Ny,Nf = size(b)
    y = zeros(Nx,Ny)
    for i in Nf:-1:1
        fistaprox!(y,psf[:,:,i],b[:,:,i],g;tol=tol, maxit=maxit,lam=lam)
        x[:,:,i] = y
    end
    return x
end

# Extent ProximalOperators to include L1_pos
import ProximalOperators: prox!

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
