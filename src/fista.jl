function fista{T}(psf::Array{T,2},b::Array{T,2},X0::Array{T,2},maxit::Int64)
    X = copy(X0)
    Xprev = X
    Y = X
    t = 1
    N = size(X,1)
    lambda = 1e-5
    psft = transpose(psf)
    obj = zeros(maxit)

    s = rand(size(psf))
    for k = 1:10
        s = imfilter(s,centered(psf))/vecnorm(s,2)
    end
    L = vecnorm(s,2)^2
    r = imfilter(Y,centered(psf),Algorithm.FFT())-b
    gradY = imfilter(r,centered(psft),Algorithm.FFT())
    k = 0
    while k < maxit
        k += 1
        obj[k] = 0.5*vecnorm(r,2)^2
        X = max.(0.0,real(Y-(1/L)*gradY))
        tnew = (1+sqrt.(1+4*t^2))/2
        Y = X + ((t-1)/tnew)*(X-Xprev)
        r = imfilter(Y,centered(psf),Algorithm.FFT())-b
        gradY = imfilter(r,centered(psft),Algorithm.FFT())
        Xprev = X
        t = tnew
    end
    return X
end

function fista{T}(psf::Array{T,3},b::Array{T,3},X0::Array{T,2},maxit::Int64)
    Y = similar(b)
    Nx,Ny,Nf = size(b)
    Xprev = copy(X0)
    obj = zeros(maxit,Nf)
    Threads.@threads for i in Nf:-1:1
        Y[:,:,i] = fista(psf[:,:,i],b[:,:,i],Xprev,maxit::Int64)
        Xprev = reshape(Y[:,:,i],Nx,Ny)
    end
    return Y
end
