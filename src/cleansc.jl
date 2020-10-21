"""
    cleanSC(E[;maxiter=50,ϕ=0.5])

CLEAN-SC algorithm for source identification and quantification optionally setting maximum iterations `maxiter` (default 50), and loop gain `ϕ` (default 0.5). 

References:
-	Sijtsma, P. (2007). CLEAN based on spatial source coherence. International Journal of Aeroacoustics, 6(4), 357–374.
"""
function cleanSC(E;maxiter=50,ϕ=0.5)
    @unpack CSM_s,steeringvec,M,N,Nf,fn = E
    x = zeros(N,Nf)
    @views @inbounds for j = 1:Nf
        _cleanSC!(x[:,j],steeringvec.arr[:,:,j],CSM_s.arr[:,:,j],maxiter,ϕ)
    end
    return FreqArray(x,E.fn)
end

function _cleanSC!(x,st,csm,maxiter,ϕ)
    sumCSM = sum(abs.(csm))
    sumDegradedCSM = sumCSM
    M,N = size(st)
    g = similar(st[:,1])
    h = similar(g)
    h_old = similar(h)
    H = similar(csm)
    dirtyMap = zeros(N)
    for _ in 1:maxiter
        AeroAcoustics.bf_col!(dirtyMap,st,csm)
        max_val,max_idx = findmax(dirtyMap)
        g .= st[:,max_idx]
        h .= g
        for _ in 1:20
            h_old .= h
            H .= h*h'
            H .= diagm(diag(H))
            h .= (csm*g/max_val + H*g) ./sqrt.(1 .+ mydot(g',H,g))
            if norm(h-h_old) .< 1e-6
                @goto cont
            end
        end
        @label cont 
        Pmax = zeros(N)
        Pmax[max_idx] = 1
        x .+= ϕ*max_val*Pmax
        csm -= ϕ*max_val*(h*h')
        csm[Diagonal(ones(Bool,M))] .= 0
        sumCSM = sum(abs.(csm))
        if sumCSM > sumDegradedCSM
            @goto stop
        end
        sumDegradedCSM = sumCSM
    end
    @label stop
    return nothing
end