import DSP.Periodograms: compute_window, forward_plan, fft2oneortwosided!, nextfastfft, arraysplit
import DSP: fftouttype

function stft!(out,s::AbstractVector{T}, n::Int=length(s)>>3, noverlap::Int=n>>1;
              onesided::Bool=eltype(s)<:Real, nfft::Int=nextfastfft(n), fs::Real=1,
              window::Union{Function,AbstractVector,Nothing}=nothing) where T
    #onesided && T <: Complex && error("cannot compute one-sided FFT of a complex signal")

    win, norm2 = compute_window(window, n)
    sig_split = arraysplit(s, n, noverlap, nfft, win)
    nout = onesided ? (nfft >> 1)+1 : nfft
    #out = zeros(stfttype(T, psdonly), nout, length(sig_split))
    tmp = Vector{fftouttype(T)}(undef, T<:Real ? (nfft >> 1)+1 : nfft)
    r = fs*norm2

    plan = forward_plan(sig_split.buf, tmp)
    offset = 0
    for sig in sig_split
        mul!(tmp, plan, sig)
        fft2oneortwosided!(out, tmp, nfft, onesided, offset)
        offset += nout
    end
    out
end

# flat_t is a helper function to reshape a S x M matrix to a vector of vectors
function flat_t(t::AbstractArray{T,N}) where {T <: AbstractFloat, N}
    if size(t, 1) < size(t, 2)
        t = permutedims(t)
    end
    return map(i -> view(t, :, i), axes(t, 2))
end

"""
    csm(t;n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n),scaling="spectrum")

Calculate cross-spectral matrix from time series `t` which is `S x M` dimensional,
where `S` is the number of samples and `M`the number of microphones.
"""
function csm(t::AbstractArray{T};n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n),scaling="spectrum",multi_thread=true) where T <: AbstractFloat
    csm(flat_t(t);n=n,noverlap=noverlap,fs=fs,win=win,scaling=scaling,multi_thread=true)
end

function csm(t::Vector{<:AbstractVector{T}};n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n),scaling="spectrum",multi_thread=true) where T <: AbstractFloat
    _foreach = AeroAcoustics.check_multithread(multi_thread)
    M = length(t)
    Nf = div(n,2)+1
    nout = div((length(t[1]) - n), n - noverlap)+1
    fc = [k*div(fs,n) for k in range(0,stop=noverlap)]
    C = Array{Complex{T}}(undef,M,M,Nf)
    S = DSP.stft.(t, n, noverlap; fs=fs, window=win, onesided=true)
    Sc = conj.(S)
    @views @inbounds _foreach(1:Nf) do ω
        for m in 1:M
            for j in m:M
                C[j, m, ω] = mean(Sc[m][ω,:] .* S[j][ω,:])
                C[m, j, ω] = conj(C[j, m, ω])  
            end
        end
    end

    if scaling == "density"
        scale = 1/(fs*sum(abs2,win))
    elseif scaling == "spectrum"
        scale = 1/(sum(win).^2)
    end

    C .*= scale
    C[:,:,2:end-1] .*= 2

    return FreqArray(C,fc)
end

"""
    csm_slow(t;n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n),scaling="spectrum")

Calculate cross-spectral matrix from time series `t` which is `S x M` dimensional,
where `S` is the number of samples and `M`the number of microphones. This version is slower than
`csm` but allocates much less (in some cases), therefore `csm_slow` could be used for long time signals.
"""
function csm_slow(t::AbstractArray{T};n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n),scaling="spectrum") where T <: AbstractFloat
    csm_slow(flat_t(t);n=n,noverlap=noverlap,fs=fs,win=win,scaling=scaling)
end

function csm_slow(t::Vector{<:AbstractVector{T}};n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n),scaling="spectrum") where T <: AbstractFloat
    M = length(t)
    Nf = div(n,2)+1
    nout = div((length(t[1]) - n), n - noverlap)+1
    weight = sum(abs2,win)
    fc = [k*div(fs,n) for k in range(0,stop=noverlap)]
    C = Array{Complex{T}}(undef,M,M,Nf)
    S = Array{Complex{T}}(undef,Nf,nout)
    Sc = similar(S)
    
    for m in 1:M
        stft!(Sc, t[m], n, noverlap; fs=fs, window=win, onesided=true)  
        conj!(Sc) 
    
        for j in m:M
            stft!(S, t[j], n, noverlap; fs=fs, window=win, onesided=true) 
            product = Sc .* S
            mean_product = mean(product, dims=2)
            C[j, m, :] .= dropdims(mean_product; dims=2)  # Drop the singleton dimension
        end
    end
    

    if scaling == "density"
        scale = 1/(fs*sum(abs2,win))
    elseif scaling == "spectrum"
        scale = 1/(sum(win).^2)
    end

    C .*= scale
    C[:,:,2:end-1] .*= 2

    for ω in 1:Nf
        C[:,:,ω] = Hermitian(C[:,:,ω],:L)
    end
    return FreqArray(C,fc)
end
