function pointspreadfunction(E::Environment{T},C::Constants,V::SteeringMatrix{T,S},cent=round.(Int,E.N/2)+1) where {T <: AbstractFloat, S <: SteeringVectorType}
    PSF = Array{T}(E.N,E.Nf)
    gm = Array{Complex{T}}(E.M)
    psf!(PSF,E,C,V,gm,cent,V.kind)
    return reshape(0.5*PSF,E.Nx,E.Ny,E.Nf)
end

function psf!(PSF,E,C,V,gm,cent,kind::Type2)
    for j in 1:length(E.f)
        gm = E.M*V.v[:,cent,j].*(E.D0[cent]./E.D[cent,:])
        for i in eachindex(E.D0)
            PSF[i,j] = abs2(V.v[:,i,j]'*gm)
        end
    end
    return PSF
end

function psf!(PSF,E,C,V,gm,cent,kind::Type3)
    for j in 1:length(E.f)
        gm = E.M*V.v[:,cent,j]
        for i in eachindex(E.D0)
            PSF[i,j] = abs2(V.v[:,i,j]'*gm)
        end
    end
    return PSF
end
