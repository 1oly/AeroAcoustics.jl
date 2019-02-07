"""
    psf(Environment[,cent])

Calculate frequency-domain point spread function using the Environment struct to access
steeringvectors. Optionally, supply the index where the psf is centered, default is (N/2)+1.
"""
function psf(E::Environment,cent::Int64=floor(Int,E.N/2)+1)
    @unpack steeringvec,M,N,Nf = E
    psf = Array{Float64, 2}(undef, N, Nf)
    for j in 1:Nf
        vc = M*view(steeringvec.arr,:,cent,j)
        for i in 1:N
            vd = @view steeringvec.arr[:,i,j]
            psf[i,j] = abs2(vd'*vc)
        end
    end
    return 0.5*psf
end
