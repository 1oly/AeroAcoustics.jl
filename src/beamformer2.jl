function beamformer2(E::Environment{T},C::Constants,V::SteeringMatrix{T}) where T <: AbstractFloat
    b = Array{T}(E.N,E.Nf)
    vd = Array{Complex{T}}(E.M)
    CSM = Array{Complex{T}}(E.M,E.M)
    # Compute beamforming
    for j in 1:length(E.f)
        CSM = E.CSM[j,:,:]
        for i in eachindex(E.D0)
            vd = V.v[:,i,j]
            b[i,j] = real(vd'*CSM*vd)
        end
    end
    return reshape(b,E.Nx,E.Ny,E.Nf)
end

function pointspreadfunction(E::Environment{T},C::Constants,V::SteeringMatrix{T}) where T <: AbstractFloat
    PSF = Array{T}(E.N,E.Nf)
    gm = Array{Complex{T}}(E.M)
    cent::Int64 = round(Int64,E.N/2)+1

    if V.kind == "II"
        for j in 1:length(E.f)
            gm = E.M*V.v[:,cent,j].*(E.D0[cent]./E.D[cent,:]) # Compensate to get correct level: Type II
            for i in eachindex(E.D0)
                PSF[i,j] = abs2(V.v[:,i,j]'*gm)
            end
        end
    elseif V.kind == "III"
        for j in 1:length(E.f)
            gm = E.M*V.v[:,cent,j]
            for i in eachindex(E.D0)
                PSF[i,j] = abs2(V.v[:,i,j]'*gm)
            end
        end
    else
        warn("Unknown steering vector type...")
        for j in 1:length(E.f)
            gm = V.v[:,cent,j]
            for i in eachindex(E.D0)
                PSF[i,j] = abs2(V.v[:,i,j]'*gm)
            end
        end
    end

    return reshape(PSF,E.Nx,E.Ny,E.Nf)
end
