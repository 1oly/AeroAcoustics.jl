function shear!(fvec,xi,xm,Ma,h)
    a = sqrt(xi[1]^2+(1-Ma^2)*(xi[2]^2+xi[3]^2))
    b = norm(xm-xi)
    fvec[1] = (1/a)*xi[1]-(1/b)*(1-Ma^2)*(xm[1]-xi[1]) - Ma
    fvec[2] = (1/a)*xi[2]-(1/b)*(xm[2]-xi[2])
    fvec[3] = xi[3] - h
end

function propagation_time(E::Environment)
    @unpack N,M,micgeom,Rxy,Ma,h,c = E
    ta = Array{Float64,2}(undef,N,M)
    xm = zeros(3,M)
    for n = 1:N
        xm[1:2,:] .= micgeom[1:2,:] .- Rxy[1:2,n]
        for m = 1:M
            res = nlsolve((fvec,x)->shear!(fvec,x,xm[:,m],Ma,h),[.1,.1,h])
            xi = res.zero
            r1 = norm(xi)
            d = xi[1]*Ma*c/r1
            c1 = d + sqrt(d^2+c^2-(Ma*c)^2)
            r2 = norm(xm[:,m]-xi)
            ta[n,m] = r1/c1 + r2/c
        end
    end
    return ta
end

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
