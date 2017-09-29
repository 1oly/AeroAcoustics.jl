function fista{T}(psf::Array{T,2},b::Array{T,2},X0::Array{T,2},maxit::Int64)
    X = copy(X0)
    Xprev = X
    Y = X
    t = 1.
    N = size(X,1)
    obj = zeros(maxit)
    psft = transpose(psf)
    Fps = fft(fftshift(psf))
    FpsT = fft(fftshift(psf'))
    s = rand(size(psf))
    for k = 1:10
        s = ifft(fft(s).*Fps)/vecnorm(s)
        #s = imfilter(s,centered(psf))/vecnorm(s,2)
    end
    L = vecnorm(s)^2
    r = real(ifft(fft(Y).*Fps))-b
    #r = imfilter(Y,centered(psf),Fill(zero(eltype(Y))),Algorithm.FFT())-b
    gradY = real(ifft(FpsT.*fft(r)))
    #gradY = imfilter(centered(psft),r,Fill(zero(eltype(r))),Algorithm.FFT())
    k = 0
    while k < maxit
        k += 1
        obj[k] = 0.5*vecnorm(r)^2
        X = max.(0.0,real(Y-(1./L)*gradY))
        tnew = 0.5*(1.+sqrt.(1.+4.*t^2))
        Y = X + ((t-1.)/tnew)*(X-Xprev)
        r = real(ifft(fft(Y).*Fps))-b
        #r = imfilter(Y,centered(psf),Fill(zero(eltype(Y))),Algorithm.FFT())-b
        gradY = real(ifft(FpsT.*fft(r)))
        #gradY = imfilter(centered(psft),r,Fill(zero(eltype(r))),Algorithm.FFT())
        Xprev = X
        t = tnew
    end
    return X,obj
end

function fista{T}(psf::Array{T,3},b::Array{T,3},X0::Array{T,2},maxit::Int64)
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

function fistalasso{T}(psf::Array{T,2},b::Array{T,2},X0::Array{T,2},maxit::Int64)
    obj = zeros(maxit)
    beta = maximum(abs.(fft(psf)))^2
    fpsf = fft(fftshift(psf))
    gamma = 1/beta
    lambda = 1e-7
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

function fistalasso{T}(psf::Array{T,3},b::Array{T,3},X0::Array{T,2},maxit::Int64)
    Y = similar(b)
    Nx,Ny,Nf = size(b)
    Xprev = copy(X0)
    obj = zeros(maxit,Nf)
    for i in 1:Nf
        Y[:,:,i],obj[:,i] = fistalasso(psf[:,:,i],b[:,:,i],Xprev,maxit)
    end

    # Warm-start:
    #Threads.@threads for i in Nf:-1:1
    #    Y[:,:,i] = fista(psf[:,:,i],b[:,:,i],Xprev,maxit)
    #    Xprev = reshape(Y[:,:,i],Nx,Ny)
    #end
    return Y,obj
end
