


f(R) =var(mean.(R))/mean(var.(R))
function f_pairs(R)
    n = length(R)
    F = zeros(n,n)
    for j in 1:n-1
        for i in j+1:n
            F[i,j] = f([R[i],R[j]])
        end
    end
    return F
end

function betadisper(X,factors, metric = Euclidean, median = true, n_perm = 1000)
    X = metric == Euclidean ? X : transform(fit(MDS,pairwise(metric(),X,X,dims = 1), distances = true))'
    levels = unique(factors)
    n = length(levels)
    centroids = []
    residuals = []
    for level in levels
        x = X[factors .== level,:]
        c = meadian == true ? DirectionalStatistics.geometric_median(x) : mean(x,dims = 1)
        r = pairwise(Euclidean(),x, c,dims = 1) .^2
        push!(centroids,c)
        push!(residuals,r)
    end 
    F = f(residuals)
    F_pairs = f_pairs(residuals)
    r = vcat(residuals...)
    perm = Vector{Float64}(undef,n_perm)
    perm_pairs = Array{Float64}(undef,n,n,n_perm)
    inds = [factors .== level for level in levels]
    for p in 1:n_perm
        shuffle!(r)
        rs = [view(r,inds[i]) for i in 1:n]
        perm[p] =f(rs)
        perm_pairs[:,:,p] = f_pairs(rs)
    end
    P = sum(perm .> F)/n_perm
    P_pairs = sum(perm_pairs .> F_pairs, dims = 3)/n_perm

    return P,P_pairs
end
