function cmf2{T,F,N}(GhG::Array{T,N},X::Array{F,N},GhPG::Array{T,N},t::Float64,lambda::Float64,kmax::Int64)
    Xprev = copy(X)
    Y = similar(X)
    k = 1
    while k<=kmax
        d = (k-1)/(k+2)
        Y .= X + d*(X-Xprev)
        Xprev .= X
        #Y = X + ((k-1)/(k+2))*(X-Xprev)
        gradY = grad(GhG,Y,GhPG,lambda)
        #@code_warntype grad(GhG,Y,GhPG,lambda)
        X = prox(Y-t*gradY)
        k += 1
    end
    return X
end

#function grad(GhG,X,GhPG,lambda)
function grad{T,N}(GhG::Array{T,N},X::Array{T,N},GhPG::Array{T,N},lambda::Float64)
    if isdiag(X)
        Y = similar(GhG)
        return Y .= diagm(-diag(GhPG-GhG*X*GhG))+lambda*eye(size(X,1))
    else
        Y = similar(GhG)
        return Y = -(GhPG-GhG*X*GhG)+lambda*eye(size(X,1))
    end
end

function prox(X)
    if isdiag(X)
        Y = similar(X)
        return Y = diagm(max.(diag(real(X)),0.0))
    else
        n = size(X,1)
        Z = similar(X)
        F = eigfact(X)
        for i in eachindex(F.values)
            F.values[i] = max.(0.0,real(F.values[i]))
        end
        for i = 1:n
            for j = 1:n
                Z[i,j] = 0.0
                for k = 1:n
                    Z[i,j] += F.vectors[i,k]*F.values[k]*F.vectors[j,k]
                end
            end
        end
        return Z
    end
end
