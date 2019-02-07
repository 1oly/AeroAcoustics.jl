"""
    steeringvectors(E)

Pre-compute steeringvectors for beamforming using an `Environment` with the needed parameters.
"""
function steeringvectors(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0 = E
    kw = 2pi*fn/c
    vi = Array{Complex{Float64}}(undef,N,M,Nf)
    for j in 1:Nf
        vi[:,:,j] .= 1 ./(D.*D0.*sum(1 ./D.^2,dims=2)).*exp.(-im.*kw[j].*(D.-D0))
    end
    E.steeringvec = FreqArray(permutedims(vi,(2,1,3)),fn)
    return E
end

function steeringvectors!(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0 = E
    kw = 2pi*fn/c
    vi = Array{Complex{Float64}}(undef,N,M,Nf)
    for j in 1:Nf
        vi[:,:,j] .= 1 ./(D.*D0.*sum(1 ./D.^2,dims=2)).*exp.(-im.*kw[j].*(D.-D0))
    end
    E.steeringvec = FreqArray(permutedims(vi,(2,1,3)),fn)
    return nothing
end
