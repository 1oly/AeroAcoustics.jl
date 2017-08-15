function fista{T}(psf::T,b::T,X0::T,maxit::Int64)
    X = copy(X0)
    Xprev = X
    Y = X
    t = 1
    N = size(X,1)
    lambda = 1e-5
    psft = transpose(psf)
    obj = zeros(maxit)
    #Fps = fft(psf)
    #FpsT = fft(transpose(psf))

    s = rand(size(psf))
    for k = 1:10
        s = imfilter(s,centered(psf))/vecnorm(s,2)
        #ifft(fft(s)*Fps)/vecnorm(s,2)
    end
    L = vecnorm(s,2)^2
    r = imfilter(Y,centered(psf))-b
    gradY = imfilter(r,centered(psft))
    k = 0
    while k < maxit
        k += 1
        obj[k] = 0.5*vecnorm(r,2)^2
        X = max.(0.0,real(Y-(1/L)*gradY))
        tnew = (1+sqrt.(1+4*t^2))/2
        Y = X + ((t-1)/tnew)*(X-Xprev)
        r = imfilter(Y,centered(psf))-b
        gradY = imfilter(r,centered(psft))
        Xprev = X
        t = tnew
    end
    return X,obj
end
