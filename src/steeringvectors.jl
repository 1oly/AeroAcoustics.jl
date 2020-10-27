ξ(Ma,a,b) = sqrt((1-Ma*a*b)^2-a^2)

function f!(F,x,Δxm,Δz2,Δz3,Ma,Δy)
    if (-0.8<=x[1]<=0.8) #&& (0.2<=abs(x[2])<=1)
        Δx20 = Δz2 * x[1]*x[2]/ξ(Ma,x[1],x[2])
        Δx21 = Δz2 * Ma*(1-Ma*x[1]*x[2])/ξ(Ma,x[1],x[2])
        Δx3  = Δz3 * x[1]*x[2]/sqrt(1-x[1]^2)
        F[1] = Δx20 + Δx21 + Δx3 - Δxm
        F[2] = x[2] - (Δxm - Δx21)/sqrt((Δxm-Δx21)^2+Δy^2)
        return nothing
    else
        return nothing
    end
end

"""
    r,pc_pm = refraction_correction(E::Environment[,h1=1.5,h2=z0-h1])

Compute propagation time correction and amplitude correction for 2D planar jet.

References:
-	Glegg, S., & Devenport, W. (2017). Aeroacoustics of low mach number flows: Fundamentals, analysis, and measurement. Chap. 10.
"""
function refraction_correction(E,h1=1.5,h2=E.z0-h1)
    r = zeros(E.N,E.M)
    pc_pm = similar(r)
    for n = 1:E.N
        for m = 1:E.M
            Δx,Δy,Δz = E.Rxy[:,n]-E.micgeom[:,m] # vector from mic to grid point
            r0 = sqrt(Δx^2+Δy^2+Δz^2)
            a0 = sqrt(Δx^2+Δy^2)/r0
            b0 = Δx/sqrt(Δx^2+Δy^2)
            Δz2,Δz3 = h1,h2
            
            res = nlsolve((F,x)->f!(F,x,Δx,Δz2,Δz3,E.Ma,Δy),[a0,b0])
            a,b = res.zero
            
            
            Δx20 = Δz2*a*b/ξ(E.Ma,a,b)
            Δx3 = Δz3*a*b/sqrt(1-a^2)
            
            r[n,m] = Δz3/sqrt(1-a^2) + Δz2*(1-E.Ma*a*b)/ξ(E.Ma,a,b)
            
            # Amplitude correction:
            pc_pi = Δz2^2/(Δz2+Δz3)^2
            pi_pt = (sqrt(1-a^2)*(1-E.Ma*a*b)^2+ξ(E.Ma,a,b))^2/(4*ξ(E.Ma,a,b)^2)
            K(R,A,D) = sqrt(a^2*(D*sqrt(1-a^2)*(1-b^2)/ξ(E.Ma,a,b) + A*(b-D/R) )^2 + a^2*(1-b^2)*(A-D*sqrt(1-a^2)*b/ξ(E.Ma,a,b))^2 + (1-a^2)*(R-D*b)^2)
            
            d = 1-E.Ma*a*b
            r2 = (1-E.Ma*a*b)*Δz2/ξ(E.Ma,a,b)
            r3 = Δz3/sqrt(1-a^2)
            ρ2 = r2/d^2
            R2 = ρ2
            R3 = r3+ρ2
            A2 = sqrt(1-a^2)*ρ2/ξ(E.Ma,a,b)
            A3 = r3+sqrt(1-a^2)*ρ2/ξ(E.Ma,a,b)
            D2 = ρ2*E.Ma*a
            D3 = ρ2*E.Ma*a
            pt_pm = R3*K(R3,A3,D3)/(R2*K(R2,A2,D2))
            pc_pm[n,m] = pc_pi * pi_pt * pt_pm
        end
    end
    return r,pc_pm
end

"""
    steeringvectors(E::Environment)

Pre-compute steeringvectors for beamforming using an `Environment` with the needed parameters.
"""
function steeringvectors(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0 = E
    kw = 2pi*fn/c
    vi = Array{ComplexF64,3}(undef,N,M,Nf)
    for j in 1:Nf
        vi[:,:,j] .= 1 ./(D.*D0.*sum(1 ./D.^2,dims=2)).*exp.(-im.*kw[j].*(D.-D0))
    end
    steeringvec = FreqArray(permutedims(vi,(2,1,3)),fn)
end

function steeringvectors!(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0,Ma,h = E
    vi = Array{ComplexF64,3}(undef,N,M,Nf)

    if E.shear
        w = 2pi*fn
        dx,pc_pm = refraction_correction(E,h)
        ta = dx./c
        for j in 1:Nf
            vi[:,:,j] .= 1 ./(ta.*c.*D0.*sum(1 ./(ta.*c).^2,dims=2)).*exp.(-im.*w[j].*ta)
        end
        if E.ampcorr
            for j in 1:Nf
                vi[:,:,j] .*= pc_pm
            end
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
