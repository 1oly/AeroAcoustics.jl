function cmf{T,C}(GhG::Array{C,2},X::Array{T,2},GhPG::Array{C,2},t::T,lambda::T,kmax::Int64)
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

function cmf{T,C}(GhG::Array{C,3},X::Array{T,2},GhPG::Array{C,3},t::T,lambda::T,kmax::Int64)
    Y = Similar(GhG)
    Nx,Ny,Nf = size(GhG)
    Xprev = copy(X)
    Threads.@threads for i in Nf:-1:1
        Y[:,:,i] = cmf(GhG[:,:,i],Xprev,GhPG[:,:,i],t,lambda,kmax)
        Xprev = reshape(Y[:,:,i],Nx,Ny)
    end
    return Y
end


function grad{T,C}(GhG::Array{C,2},X::Array{T,2},GhPG::Array{C,2},lambda::T)
    if isdiag(X)
        Y = similar(GhG)
        return Y .= diagm(-diag(GhPG-GhG*X*GhG))+lambda*eye(size(X,1))
    else
        Y = similar(GhG)
        return Y = -(GhPG-GhG*X*GhG)+lambda*eye(size(X,1))
    end
end

function prox{T}(X::Array{T,2})
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
