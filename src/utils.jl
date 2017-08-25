function SPL{T}(p::Array{T})
    s = similar(p)
    s[p.>0] = 10*log10.(p[p.>0]/4e-10)
    s[p.<0] = -350
    return s
end
SPL(p::Number) = 10*log10(p/4e-10)

function shear(xi,fvec,xm,M,h)
    a = sqrt(xi[1]^2+(1-M^2)*(xi[2]^2+xi[3]^2))
    b = sqrt((xm[1]-xi[1])^2+(xm[2]-xi[2])^2+(xm[3]-xi[3])^2)
    fvec[1] = (1/a)*xi[1]-(1/b)*(1-M^2)*(xm[1]-xi[1]) - M
    fvec[2] = (1/a)*xi[2]-(1/b)*(xm[2]-xi[2])
    fvec[3] = xi[3] - h
end
