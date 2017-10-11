function fista(psf::Array{T,2},b::Array{T,2},X0::Array{T,2},maxit::Int64) where T
    X = copy(X0)
    Xprev = X
    Y = X
    t = 1.
    Nx,Ny = size(X)
    obj = zeros(maxit)
    psft = transpose(psf)
    center = round(Int64,Nx/2)+1,round(Int64,Ny/2)+1
    Fps = fft(circshift(psf,center))
    FpsT = fft(circshift(psf',center))
    s = rand(size(psf))
    for k = 1:10
        s = ifft(fft(s).*Fps)/vecnorm(s)
        #s = imfilter(s,centered(psf))/vecnorm(s)
    end
    L = vecnorm(s)^2
    r = real(ifft(Fps.*fft(Y)))-b
    #r = imfilter(Y,centered(psf),Fill(zero(eltype(Y))),Algorithm.FFT())-b
    gradY = real(ifft(FpsT.*fft(r)))
    #gradY = imfilter(r,centered(psft),Fill(zero(eltype(r))),Algorithm.FFT())
    k = 0
    while k < maxit
        k += 1
        obj[k] = 0.5*vecnorm(r)^2
        X = max.(0.0,real(Y-(1./L)*gradY))
        tnew = 0.5*(1.+sqrt.(1.+4.*t^2))
        Y = X + ((t-1.)/tnew)*(X-Xprev)
        r = real(ifft(Fps.*fft(Y)))-b
        #r = imfilter(Y,centered(psf),Fill(zero(eltype(Y))),Algorithm.FFT())-b
        gradY = real(ifft(FpsT.*fft(r)))
        #gradY = imfilter(r,centered(psft),Fill(zero(eltype(r))),Algorithm.FFT())
        Xprev = X
        t = tnew
    end
    return X,obj
end

function fista(psf::Array{T,3},b::Array{T,3},X0::Array{T,2},maxit::Int64) where T
    Y = similar(b)
    Nx,Ny,Nf = size(b)
    Xprev = copy(X0)
    obj = zeros(maxit,Nf)
    for i in 1:Nf
        Y[:,:,i],obj[:,i] = fista(psf[:,:,i],b[:,:,i],Xprev,maxit)
    end

    # Warm-start:
    #Threads.@threads for i in Nf:-1:1
    #    Y[:,:,i] = fista(psf[:,:,i],b[:,:,i],Xprev,maxit)
    #    Xprev = reshape(Y[:,:,i],Nx,Ny)
    #end
    return Y,obj
end

function fistalasso(psf::Array{T,2},b::Array{T,2},X0::Array{T,2},maxit::Int64,lambda::T) where T
    obj = zeros(maxit)
    Nx,Ny = size(X0)
    beta = maximum(abs.(fft(psf)))^2
    center = round(Int64,Nx/2)+1,round(Int64,Ny/2)+1
    fpsf = fft(circshift(psf,center))
    #fpsf = fft(fftshift(psf))
    gamma = 1/beta
    prox(x,gamma,lambda) = x - x./max.(abs.(x)./(lambda.*gamma), 1)
    a = 10
    X = b
    Z = X
    for k = 1:maxit
        Xprev = X
        gradf = real(ifft(fpsf'.*fft(ifft(fpsf.*fft(Z))-b)))
        X = prox(Z-gamma.*gradf,gamma,lambda)
        #X = max.(0.0,Z-gamma.*gradf)
        alpha = k/(k+1+a)
        Z = X + alpha.*(X-Xprev)
        obj[k] = 0.5*vecnorm(ifft(fpsf.*fft(Z))-b)^2
    end
    return X,obj
end

function fistalasso(psf::Array{T,3},b::Array{T,3},X0::Array{T,2},maxit::Int64,lambda::T) where T
    Y = similar(b)
    Nx,Ny,Nf = size(b)
    Xprev = copy(X0)
    obj = zeros(maxit,Nf)
    for i in 1:Nf
        Y[:,:,i],obj[:,i] = fistalasso(psf[:,:,i],b[:,:,i],Xprev,maxit,lambda)
    end

    # Warm-start:
    #Threads.@threads for i in Nf:-1:1
    #    Y[:,:,i] = fista(psf[:,:,i],b[:,:,i],Xprev,maxit)
    #    Xprev = reshape(Y[:,:,i],Nx,Ny)
    #end
    return Y,obj
end
