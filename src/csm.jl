function csm(t::AbstractArray{T}) where T <: AbstractFloat
    const n = 1024
    const noverlap = div(n,2)
    const fs = 51200
    const win = DSP.hanning(n)
    const nout = div((size(t,1) - n), n - noverlap)+1
    const M = size(t,2)
    const Nf = div(n,2)+1
    const weight = dot(win,win)

    s = zeros(Complex{T},Nf,nout,M)
    C = zeros(Complex{T},Nf,M,M)

    for i in 1:size(t,2)
        s[:,:,i] = DSP.stft(t[:,i], n, noverlap; fs=fs, window=win, onesided=true)
    end

    for l = 1:nout
        for ω = 1:Nf
            C[ω,:,:] += s[ω,l,:]*s[ω,l,:]'
        end
    end
    C .*= 2/n/weight/nout;
    return C
end
