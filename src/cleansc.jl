"""
    cleanSC(E[;maxiter=50,ϕ=0.5,stopn=10,peak_removal=false,trust_region=nothing])

CLEAN-SC algorithm for source identification and quantification optionally setting maximum iterations `maxiter` (default 50), and loop gain `ϕ` (default 0.5). Additionally, a stopping criterion `max_peak[i] > max_peak[i-stopn]` can be set with parameter `stopn`. Subtraction of peak sources in the dirty map is supported by `peak_removal=true/false` to remove peak sources outside `trust_region` with limits given as `[xmin,xmax,ymin,ymax]` e.g. `trust_region=[-0.75;1.5;-1.0;1.0]`.

References:
-	Sijtsma, P. (2007). CLEAN based on spatial source coherence. International Journal of Aeroacoustics, 6(4), 357–374.

With inspiration from https://github.com/acoular/acoular/blob/66cba3cffb3bc72602c869f99347be76798f4ac1/acoular/fbeamform.py#L1496
"""
function cleanSC(E;maxiter=50,ϕ=0.5,stopn=10,peak_removal=false,trust_region=nothing)
    @unpack CSM_s,steeringvec,M,N,Nf,fn = E
    x = zeros(N,Nf)
    
    if peak_removal
        xi = findall((E.rx.>=trust_region[1]) .& (E.rx.<=trust_region[2]))
        yi = findall((E.ry.>=trust_region[3]) .& (E.ry.<=trust_region[4]))
        trust_indices = LinearIndices((E.Nx,E.Ny))[xi,yi]
    else
        trust_indices = LinearIndices((E.Nx,E.Ny))
    end
    
    @views @inbounds for j = 1:Nf
        _cleanSC!(x[:,j],steeringvec.arr[:,:,j],CSM_s.arr[:,:,j],maxiter,ϕ,stopn,peak_removal,trust_indices)
    end
    return FreqArray(x,E.fn)
end

function _cleanSC!(x,st,csm,maxiter,ϕ,stopn,peak_removal,trust_indices)
    M,N = size(st)
    g = similar(st[:,1])
    h = similar(g)
    h_old = similar(h)
    H = similar(csm)
    max_val = zeros(maxiter)
    dirtyMap = zeros(N)
    for i in 1:maxiter
        AeroAcoustics.bf_col!(dirtyMap,st,csm)
        # Numerical issue with Amiet correction can cause NaNs
        if isempty(filter(!isnan,dirtyMap)) 
            @goto stop
        else
            max_val[i],max_idx = findmax(filter(!isnan,dirtyMap))
        end
        g .= st[:,max_idx]
        h .= g
        for _ in 1:20
            h_old .= h
            H .= h*h'
            H .= diagm(diag(H))
            h .= (csm*g/max_val[i] + H*g) ./sqrt.(1 .+ mydot(g',H,g))
        end
        Pmax = zeros(N)
        Pmax[max_idx] = 1
        
        if peak_removal && max_idx ∈ trust_indices
            x .+= ϕ*max_val[i]*Pmax
            csm -= ϕ*max_val[i]*(h*h')
        else
            csm -= max_val[i]*(h*h')
        end
        
        csm[Diagonal(ones(Bool,M))] .= 0
        if (i > stopn) && (max_val[i] > max_val[i-stopn])
            @goto stop
        end
    end
    @label stop
    return nothing
end