# structs

struct Disp
    F ::Float64
    pairwise_F ::Array{Float64} 
    medians ::Union{Tuple,Vector}
    residuals ::Vector{Matrix{Float64}} 
    means ::Vector{Float64}
    group  ::Vector
    levels  ::Vector
end
#Helper functions

# calculate distances from medians
function Resids(x,c)
    d = x .-c
    sum(d .^2,dims =2) 
end

# doube center matrix
function dblcen(D)
    A = D .- mean(D,dims = 1)
    return A .- mean(A,dims = 2) 
end

# calculate medians 
# Get spatial median as a row vector
spatial_median(X ::Array) = DirectionalStatistics.geometric_median([row for row in eachrow(X)])'
#apply `spatial_median` function to get pos and neg median for each group
function spatialMed(vectors, group, pos)

    spMedPos = [spatial_median(vectors[group .== g,pos]) for g in unique(group)]
    spMedNeg = [spatial_median(vectors[group .== g, .!pos]) for g in unique(group)]

    return (spMedPos,spMedNeg )
end
function spatialMed(vectors, group)

    return [spatial_median(vectors[group .== g,:]) for g in unique(group)]
end

function dispersion(D,group;metric = false)
    """
    Quantify the dispersion around the group spatial median for each group in D, according to the given distance metric.
    See Distances.jl for list of available metrics and the process of implementic your own metric. Once X is converted 
    to a distance matrix, it is then transformed by principal coordinate analysis before medians and distances are calculated.
    The algorithm is based on Anderson (2006). But we aim to be consistent with the betadisper implementation in the
    R package Vegan, which has implemented some changes since the original paper was published.
    
    The function accepts:
        D: A dissimiarity or distance matrix
        group: A vector containing group labels
        metric: (key word) boolean value specifying whether or not the distance used used metric
    
    This function returns a named tuple containing:
        F = Global F-Statistic 
        pairwise_F = pairwise F-statistics
        medians = spatial (aka geometric) median of each group in the transformed coordinates
        residuals = vector containing each observations distance from its group median
        means = vector containing the means distance to median for each group
        group = the original grouping vector
        levels = unique items from 'group'
    
    Anderson, M.J. (2006) Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62(1), 245--253.
    https://github.com/vegandevs/vegan/blob/master/R/betadisper.R
    """
    levels = sort(unique(group))
    # vegdist objects contain lower triangular matrices wheras D is symmetric
    # thus we skip a D + D' step here.
    if metric == false
        # We double center the matrix as in vegan.betadisper. Though this is not identical Anderson (2006).
        A = dblcen(D .^2 ) 
        e = eigen(Hermitian(-A/2))
        vectors = e.vectors
        eig = e.values 
        # Anderson (2006) multiplies eigen vectors corresponding to negative eigen values by sqrt(-1).
        # We ommit this step in keeping with vegan.betadisper
        vectors = vectors * diagm(sqrt.(abs.(eig))) 
        # record indices corresponding to postitive eigenvalues
        pos = real.(eig) .> 0.0 
        medians = spatialMed(vectors, group, pos)
        # calculate distances (to median) for "pos" and "neg" vectors separately. in orer to 
        # subtract "neg" from "pos". See Anderson (2006) for details.
        dis_pos = [Resids(vectors[group .== levels[i],pos], medians[1][i]) for i in 1:length(levels)]
        dis_neg = [Resids(vectors[group .== levels[i], .!pos], medians[2][i]) for i in 1:length(levels)]
        # Where dis_neg > dis_pos, set residual = 0 by taking only the real part 
        # of the square root of a negative. This was implemented in vegan.betadisper after discussion in
        # issue #306. But this is not done in Anderson (2006)
        residuals = [(real.(sqrt.(Complex.(dis_pos[i] .- dis_neg[i]))) )   for i in 1:length(dis_pos)] 
    else
        medians = spatialMed(D, group)
        residuals = [Resids(D[group .== levels[i],:], medians[i]) for i in 1:length(levels)]
    end

    F = f(residuals)
    F_pairs = f_pairs(residuals)
    means = NamedArray(mean.(residuals),string.(levels),"group")
    return Disp(F, F_pairs, medians, residuals, means,group, levels)
end


function dispersion(X,group, distance)
    """
        Quantify the dispersion around the group spatial median for each group in X, according to the given distance metric.
        See Distances.jl for list of available metrics and the process of implementic your own metric. Once X is converted 
        to a distance matrix, it is then transformed by principal coordinate analysis before medians and distances are calculated.
        The algorithm is based on Anderson (2006). But we aim to be consistent with the betadisper implementation in the
        R package Vegan, which has implemented some changes since the original paper was published.
        
        The function accepts:
            X: A numerical array
            group: A vector containing group labels
            distance: A distance/dissimilarity measure implemented in Distances.jl
        
            This function returns a named tuple containing:
            F = Global F-Statistic 
            pairwise_F = pairwise F-statistics
            medians = spatial (aka geometric) median of each group in the transformed coordinates
            residuals = vector containing each observations distance from its group median
            means = vector containing the means distance to median for each group
            group = the original grouping vector
            levels = unique items from 'group'
        
        Anderson, M.J. (2006) Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62(1), 245--253.
        https://github.com/vegandevs/vegan/blob/master/R/betadisper.R
    """

    D = Distances.pairwise(distance(),X,X,dims = 1)
    return dispersion(D,group)
  

end