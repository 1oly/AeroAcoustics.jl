function fista(psf::Array{T,2},b::Array{T,2},X0::Array{T,2};maxit::Int64=1000,tol::T=1e-6) where T <: AbstractFloat
    X = copy(X0)
    Xprev = X
    const Nx,Ny = size(X)
    obj = zeros(maxit)
    center = round(Int64,Nx/2)+1,round(Int64,Ny/2)+1
    Fps = fft(circshift(psf,center))
    FpsT = fft(circshift(psf',center))
    s = rand(size(psf))
    for k = 1:10        # TODO: Try powm from IterativeSolvers.jl
        s = ifft(fft(s).*Fps)/vecnorm(s)
    end
    const L = vecnorm(s)^2
    const beta = 1 ./L

    it = 0
    while it < maxit
        it += 1
        x_extr = X + (it-2)/(it+1)*(X - Xprev)
        res = real(ifft(Fps.*fft(x_extr)))-b
        y = real(x_extr-beta*real(ifft(FpsT.*fft(res))))
        Xprev .= X
        X = nonneg!(y)
        if norm(x_extr-X, Inf)/beta <= tol*(1+norm(X, Inf))
            break
        end
    end
    return X, it
end

function fista(psf::Array{T,3},b::Array{T,3},X0::Array{T,2};maxit::Int64=1000,tol::T=1e-8,warmstart::Bool=true) where T <: AbstractFloat
    Y = similar(b)
    const Nx,Ny,Nf = size(b)
    Xprev = copy(X0)
    #obj = zeros(maxit,Nf)
    d = DeconvResults(Symbol("fista"))
    if warmstart
        Ind = Nf:-1:1
    else
        Ind = 1:Nf
    end

    for i in Ind
        Y[:,:,i],it = fista(psf[:,:,i],b[:,:,i],Xprev;maxit=maxit,tol=tol)
        warmstart ? Xprev = reshape(Y[:,:,i],Nx,Ny) : nothing
        push!(d.maxit,it)
    end
    d.xopt, d.tol = Y, tol
    d.maxit = warmstart ? reverse!(d.maxit) : d.maxit
    return d
end
