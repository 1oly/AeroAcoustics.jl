"""
    fista(env::Environment,b ,p [,f<:AbstractArray]; tol=1e-6, maxit=1000)
Performs deconvolution with FISTA with a nonnegativity constraint for all frequency bins in `env.fn` or frequencies in `f`. `f` must contain values that match exact a (sub)set of the values in `env.fn`. 
Refs: 
    Beck, A., & Teboulle, M. (2009). A fast iterative shrinkage-thresholding algorithm for linear inverse problems. 
        SIAM journal on imaging sciences, 2(1), 183-202.
    Lylloff, O., Fernández-Grande, E., Agerkvist, F., Hald, J., Tiana Roig, E., & Andersen, M. S. (2015). Improving the efficiency of deconvolution algorithms for sound source localization. 
        The journal of the acoustical society of America, 138(1), 172-180
"""

function fista(env::Environment,b::FA,p::FA,f::T=env.fn; tol::Real=1e-6, maxit::Int=1000) where {T <: AbstractArray, FA <: FreqArray}
    env.fn == f ? (f_inds=1:env.Nf) : (f_inds = findall(x->x in f, env.fn))
    x = zeros(size(b,1),length(f_inds))
    center = ceil(Int64,env.N / 2) 
    for (i,fn) in enumerate(f_inds)
        pp = power_iteration(p[:,fn])
        L = norm(pp)^2
        gam = 1/L
        @views _fista!(x[:,i], b[:,fn], circshift(p[:,fn],center),gam; maxit=maxit, tol=tol)
    end
    return FreqArray(x,env.fn[f_inds])
end

function _fista!(x,b,psf,gam;maxit=1000,tol=1e-6)
    x_prev = copy(x)
    Ax = zeros(ComplexF64,size(b))
    Fifft = plan_ifft(Ax)
    Fpsf = fft(psf)
    for it = 1:maxit
        x_ex = x + (it-2)/(it+1)*(x - x_prev)
        Ax = Fifft*(FFTW.fft(x_ex).*Fpsf)
        res = Ax - b
        Atx = Fifft*(FFTW.fft(res).*Fpsf)
        y = x_ex - gam*(Atx)
        x_prev .= x
        x .= max.(real.(y),0.0)
        # stopping criterion
        if norm(x_ex-x)/gam <= tol*(1+norm(x))
            println("breaking at iter = $(it)...")
            break
        end
    end
    return nothing
end

function power_iteration(A)
    x = rand(ComplexF64,size(A))
    for i = 1:10
        x = FFTW.ifft(FFTW.fft(x).*fft(A))/norm(x)
    end
    return x
end