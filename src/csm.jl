function csm(t::AbstractArray{T}) where T <: AbstractFloat
    const n = 1024
    const noverlap = div(n,2)
    const fs = 51200
    const win = DSP.hanning(n)
    const nout = div((size(t,1) - n), n - noverlap)+1
    const M = size(t,2)
    const Nf = div(n,2)+1
    const weight = dot(win,win)

    Pxy = zeros(Complex{T},Nf)
    C = zeros(Complex{T},Nf,M,M)

    for m in 1:M
        for j in m:M
            Pxy = mean(conj(DSP.stft(t[:,m], n, noverlap; fs=fs, window=win, onesided=true)).*DSP.stft(t[:,j], n, noverlap; fs=fs, window=win, onesided=true),2)
            for ω in 1:Nf
                C[ω,j,m] = Pxy[ω]
            end
        end
    end

    C .*= 2/n/weight

    for ω in 1:Nf
        C[ω,:,:] = Hermitian(C[ω,:,:],:L)
    end
    return C
end
