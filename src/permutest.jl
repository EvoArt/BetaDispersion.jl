
function permutest(disp,n_perm = 1000)
    """
    This function permutes the residuals returned by 'dispersion' to provide global and pairwise
    P-values. Currently, no corrections are made for multiple testing. MultipleTesting.jl provides
    this functionality, if desired

    The function accepts:
        disp: The named tuple returned by 'dispersion'
        n_perm: the number of permutations

    This function returns a named tuple containing:
        P = Global P-value
        pairwise_P = pairwise P-value
        F = Global F-Statistic 
        pairwise_F = pairwise F-statistics

    https://github.com/juliangehring/MultipleTesting.jl
    """
    # read in values from `disp`
    F = disp.F
    R = disp.residuals
    F_pairs = disp.pairwise_F
    group = disp.group
    levels = disp.levels
    level_names = string.(levels)
    inds = [group .== level for level in levels]
    n = length(levels)

    # pre-calculate values for ANOVA function
    r = vcat(R...)
    X = vcat(R...)
    N = length(X)
    nj = Tuple(length.(R))
    X̅ = mean(X)
    k = length(R)
    N_p,nj_p,X̅_p = get_pars(R)

    # generate empty containers for F values
    perm = Vector{Float64}(undef,n_perm)
    perm_pairs = Array{Float64}(undef,n,n,n_perm)
    
    # run permutation
    for p in 1:n_perm
        shuffle!(r)
        rs = [view(r,inds[i]) for i in 1:n]
        perm[p] =f(rs,N,nj,X̅,k)
        perm_pairs[:,:,p] = f_pairs(rs,N_p,nj_p,X̅_p)
    end
    # calculate P-values and return P and F values
    P = sum(perm .> F)/n_perm
    P_pairs = NamedArray(LowerTriangular(sum(perm_pairs .> F_pairs, dims = 3)[:,:,1]/n_perm),( level_names,level_names), ("group","group"))
    F_pairs = NamedArray(LowerTriangular(F_pairs),( level_names,level_names ), ("group","group"))
    
    return (P =P,pairwise_P = P_pairs,F = F, pairwise_F =F_pairs)
end


