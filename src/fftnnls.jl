function fftnnls(psf::Array{T,2},b::Array{T,2},X0::Array{T,2},maxit::Int64) where T <: AbstractFloat
    X = copy(X0)
    Nx,Ny = size(X)
    obj = zeros(maxit)
    center = round(Int64,Nx/2)+1,round(Int64,Ny/2)+1
    Fps = fft(circshift(psf,center))
    FpsT = fft(circshift(psf',center))
	alpha = 0.0
    k = 0
    while k < maxit
        k += 1
		r = real(ifft(Fps.*fft(X)))-b
		obj[k] = 0.5*sum(abs2,r) #0.5*vecnorm(r)^2
	    g = real(ifft(FpsT.*fft(r)))
		w = g
		w[(X .== 0.0) .& (w .> 0.0)] = 0.0
		g = real(ifft(FpsT.*fft(w)))
		alpha = dot(vec(g),vec(r))/dot(vec(g),vec(g))
		X = max.(0.0,real(X-alpha*w))
    end
    return X,obj
end

function fftnnls(psf::Array{T,3},b::Array{T,3},X0::Array{T,2},maxit::Int64) where T <: AbstractFloat
    X = similar(b)
    Nx,Ny,Nf = size(b)
    obj = zeros(maxit,Nf)
    for i in 1:Nf
        X[:,:,i],obj[:,i] = fftnnls(psf[:,:,i],b[:,:,i],X0,maxit)
    end

    # Warm-start:
    #Threads.@threads for i in Nf:-1:1
    #    Y[:,:,i] = fista(psf[:,:,i],b[:,:,i],Xprev,maxit)
    #    Xprev = reshape(Y[:,:,i],Nx,Ny)
    #end
    return X,obj
end
