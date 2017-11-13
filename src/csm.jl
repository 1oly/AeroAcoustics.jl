function csm(t::AbstractArray{T}) where T <: AbstractFloat
    const n = 1024
    const noverlap = div(n,2)
    const fs = 51200
    const win = DSP.hanning(n)
    const nout = div((size(t,1) - n), n - noverlap)+1
    const M = size(t,2)
    const Nf = div(n,2)+1
    const weight = dot(win,win)
    const df = div(fs,n)
    const nspec = div(n,2)

    Pxy = Array{Complex{T}}(Nf)
    C = Array{Complex{T}}(Nf,M,M)

    for m in 1:M
        for j in m:M
            mean!(Pxy,conj!(DSP.stft(t[:,m], n, noverlap; fs=fs, window=win, onesided=true)).*DSP.stft(t[:,j], n, noverlap; fs=fs, window=win, onesided=true))
            for ω in 1:Nf
                C[ω,j,m] = Pxy[ω]
            end
        end
    end

    C .*= 2/n/weight

    for ω in 1:Nf
        C[ω,:,:] = Hermitian(C[ω,:,:],:L)
    end

    return CrossSpectralMatrix(real.(C),imag.(C),Array{T,1}([k*df for k in 0:nspec]),false)
    #return C
end

# TODO: diagrm! Not finished!!
function diagrm!(csm::CrossSpectralMatrix{T}) where T <: AbstractFloat
    Nf,M,M = size(csm.csmReal)
    m = JuMP.Model(solver = SCS.SCSSolver())
    JuMP.@variable(m,d[1:2M])
    JuMP.@objective(m,Min,sum(d))
    D = diagm(d)
    for ω in 1:Nf
        # TODO: Make imag_expand!() and call that here:
        C = [csm.csmReal[ω,:,:] -csm.csmImag[ω,:,:]; csm.csmImag[ω,:,:] csm.csmReal[ω,:,:]]
        JuMP.@SDconstraint(m,(C+D>=0))
        JuMP.solve(m)
        dopt = JuMP.getvalue(d)
        for i in 1:M
            csm.csmReal[ω,i,i] += dopt[i]
        end
        #csm.csmImag[ω,:,:] += diagm(dopt[M+1:2M])
    end
    return CrossSpectralMatrix(csm.csmReal,csm.csmImag,csm.binCenterFrequenciesHz,true)
end
