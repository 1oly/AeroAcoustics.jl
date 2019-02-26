"""
    psf(Environment[,cent])

Calculate frequency-domain point spread function using the Environment struct to access
steeringvectors. Optionally, supply the index where the psf is centered, default is (N/2)+1.
"""
function psf(E::Environment,cent::Int64=floor(Int,E.N/2)+1)
    @unpack steeringvec,M,N,Nf,fn = E
    psf = Array{Float64, 2}(undef, N, Nf)
    for j in 1:Nf
        vc = M*view(steeringvec.arr,:,cent,j)
        for i in 1:N
            vd = @view steeringvec.arr[:,i,j]
            psf[i,j] = abs2(vd'*vc)
        end
    end
    return FreqArray(0.5*psf,fn)
end

"""
    psf_col!(p,steeringvec,cent)

Calculate single frequency point spread function as a column vector. Used for DAMAS.
"""
function psf_col!(p,steeringvec,cent)
    M = size(steeringvec,1)
    vc = M*view(steeringvec,:,cent)
    for i in eachindex(p)
        p[i] = 0.5*abs2(view(steeringvec,:,i)'*vc)
    end
    return p
end
