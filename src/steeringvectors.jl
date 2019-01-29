"""
    steeringvectors(E)

Pre-compute steeringvectors for beamforming using an `Environment` with the needed parameters.
"""
function steeringvectors(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0 = E
    kw = collect(2pi*fn/c)
    vi = Array{Complex{Float64}}(undef,M,N,Nf)
    Dsum = sum(1 ./D.^2,dims=2)
    for j in 1:Nf
        for i in eachindex(D0)
            for m in 1:M
                vi[m,i,j] = 1/(D[i,m]*D0[i]*Dsum[i])*exp(-im*kw[j]*(D[i,m]-D0[i]))
            end
        end
    end
    E.steeringvec = FreqArray(vi,fn)
    return E
end

function steeringvectors!(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0 = E
    kw = collect(2pi*fn/c)
    vi = Array{Complex{Float64}}(undef,M,N,Nf)
    Dsum = sum(1 ./D.^2,dims=2)
    for j in 1:Nf
        for i in eachindex(D0)
            for m in 1:M
                vi[m,i,j] = 1/(D[i,m]*D0[i]*Dsum[i])*exp(-im*kw[j]*(D[i,m]-D0[i]))
            end
        end
    end
    E.steeringvec = FreqArray(vi,fn)
    return nothing
end
