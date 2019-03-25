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
    steeringvec = FreqArray(permutedims(vi,(2,1,3)),fn)
end

function steeringvectors!(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0,Ma,h = E
    vi = Array{Complex{Float64}}(undef,N,M,Nf)

    if Ma != 0.0
        w = 2pi*fn
        ta = AeroAcoustics.propagation_time(E)
        for j in 1:Nf
            vi[:,:,j] .= (1 ./M).*(ta.*c./D0).*exp.(-im.*w[j].*ta)
        end
    else
        kw = 2pi*fn/c
        for j in 1:Nf
            vi[:,:,j] .= 1 ./(D.*D0.*sum(1 ./D.^2,dims=2)).*exp.(-im.*kw[j].*(D.-D0))
        end
    end

    E.steeringvec = FreqArray(permutedims(vi,(2,1,3)),fn)
    return nothing
end
