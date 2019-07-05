# flat_t is a helper function to reshape a S x M matrix to a vector of vectors
function flat_t(t::AbstractArray{T,N}) where {T <: AbstractFloat, N}
    if size(t,1) < size(t,2)
        t = permutedims(t)
    end
    ta = Array{Vector{T}}(undef,size(t,2))
    for i in axes(t,2)
        ta[i] = view(t,:,i)
    end
    return ta
end

"""
    csm(t;n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n))

Calculate cross-spectral matrix from time series `t` which is `S x M` dimensional,
where `S` is the number of samples and `M`the number of microphones.
"""
function csm(t::AbstractArray{T};n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n)) where T <: AbstractFloat
    csm(flat_t(t);n=n,noverlap=noverlap,fs=fs,win=win)
end

function csm(t::Vector{Vector{T}};n=1024,noverlap=div(n,2),fs=1,win=DSP.hanning(n)) where T <: AbstractFloat
    M = length(t)
    Nf = div(n,2)+1
    nout = div((length(t[1]) - n), n - noverlap)+1
    weight = sum(abs2,win)
    fc = [k*div(fs,n) for k in range(0,stop=noverlap)]
    C = Array{Complex{T}}(undef,M,M,Nf)
    S = DSP.stft.(t, n, noverlap; fs=fs, window=win, onesided=true)
    Sc = conj.(S)
    for m in 1:M
        for j in m:M
            C[j,m,:] .= dropdims(mean(LazyArray(@~ Sc[m].*S[j]),dims=2);dims=2)
        end
    end

    C .*= 2/n/weight

    for ω in 1:Nf
        C[:,:,ω] = Hermitian(C[:,:,ω],:L)
    end
    return FreqArray(C,fc)
end
