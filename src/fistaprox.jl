function fistaprox!(x::Array{T,2}, psf::Array{T,2}, b::Array{T,2}, g; tol::R=1e-6, maxit::Int=5000,lambda::R=0.01) where {T <: AbstractFloat, R <: Real}
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
    const beta = 1 ./L
    #lambda = lam*norm(real(ifft(FpsT.*fft(b))),Inf)
    f = g(lambda)
    it = 0
    while it < maxit
        it += 1
        x_extr = x + (it-2)/(it+1)*(x - x_prev)          # extrapolation step
        res = real(ifft(Fps.*fft(x_extr)))-b             # res = A*x_extr - b
        y = real(x_extr-beta*real(ifft(FpsT.*fft(res)))) # y = x_extr - gam*(A'*res)
        x_prev .= x                                      # store current iterate
        prox!(x, f, y, beta)
        nonneg!(x)
        # stopping criterion
        if norm(x_extr-x, Inf)/beta <= tol*(1+norm(x, Inf))
            break
        end
    end
    return x,it #TODO: Implement DeconvResults for 2D...
end

# For multidimensional arrays, this is per definition warm-starting:
function fistaprox!(x::Array{T,3}, psf::Array{T,3}, b::Array{T,3}, g; tol::R=1e-3, maxit::Int=5000,lambda::R=0.01,warmstart::Bool=true) where {T <: AbstractFloat, R <: Real}
    const Nx,Ny,Nf = size(b)
    y = zeros(Nx,Ny)
    d = DeconvResults(Symbol(g))
    if warmstart
        Ind = Nf:-1:1
    else
        Ind = 1:Nf
    end
    for i in Ind
        _,it = fistaprox!(y,psf[:,:,i],b[:,:,i],g;tol=tol, maxit=maxit,lambda=lambda)
        x[:,:,i] = y
        push!(d.maxit,it)
    end
    d.xopt, d.tol, d.lambda = x, tol, lambda
    d.maxit = warmstart ? reverse!(d.maxit) : d.maxit
    return d
end

function nonneg!(x::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    for k in eachindex(x)
        if x[k] < 0
            x[k] = zero(T)
        end
    end
    return x
end
