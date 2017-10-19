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

    return reshape(0.5*PSF,E.Nx,E.Ny,E.Nf)
end
