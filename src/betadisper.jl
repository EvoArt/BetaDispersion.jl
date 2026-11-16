"""
    BetaDisper

Result of [`betadisper`](@ref).

- `distances`: distance of each observation to its group's centre
- `group`, `levels`: group of each observation, and the sorted group levels
- `group_distances`: mean distance to centre in each group (in `levels` order)
- `centroids`: group centres (spatial medians or centroids), one row per level,
  in the principal coordinate space
- `vectors`, `eig`: principal coordinates of the observations and their
  eigenvalues (zero eigenvalues removed; negative ones kept, as in vegan)
- `type`: `:median` or `:centroid`
- `labels`: observation names (may be empty)
"""
struct BetaDisper{G}
    distances::Vector{Float64}
    group::Vector{G}
    levels::Vector{G}
    group_distances::Vector{Float64}
    centroids::Matrix{Float64}
    vectors::Matrix{Float64}
    eig::Vector{Float64}
    type::Symbol
    bias_adjust::Bool
    labels::Vector{String}
end

"""
    betadisper(D, group; type=:median, bias_adjust=false, sqrt_dist=false,
               add=false, labels=String[])
    betadisper(X, group, metric; kw...)

Distances of observations to their group centre in principal coordinate
space (Anderson 2006), the same as `vegan::betadisper(d, group)`.

`D` is a symmetric dissimilarity matrix; or give a data matrix `X` (one row
per observation) and a Distances.jl metric such as `BrayCurtis()`.
`type=:centroid` uses group centroids instead of spatial medians.
`add=:lingoes` or `:cailliez` corrects non-Euclidean dissimilarities;
`sqrt_dist=true` takes square roots first. `bias_adjust=true` corrects for
small group sizes. Observations whose group is `missing` are dropped.

Test the differences with [`permutest`](@ref) or [`dispersion_anova`](@ref).
"""
function betadisper(D::AbstractMatrix, group::AbstractVector; type::Symbol=:median,
                    bias_adjust::Bool=false, sqrt_dist::Bool=false, add=false, labels=String[])
    type in (:median, :centroid) || throw(ArgumentError("type must be :median or :centroid, got :$type"))
    n = size(D, 1)
    size(D, 2) == n || throw(ArgumentError("D must be a square dissimilarity matrix, got $(size(D))"))
    length(group) == n || throw(DimensionMismatch("group has $(length(group)) entries but D has $n rows"))
    labs = isempty(labels) ? String[] : string.(collect(labels))
    keep = findall(!ismissing, group)
    if length(keep) < n
        @info "$(n - length(keep)) observation(s) with missing group removed"
        D = D[keep, keep]
        labs = isempty(labs) ? labs : labs[keep]
    end
    g = [x for x in group[keep]]                  # drops Missing from the eltype
    TOL = sqrt(eps())
    any(<(-TOL), D) && throw(ArgumentError("dissimilarities must be non-negative"))
    d = Matrix{Float64}(D)
    sqrt_dist && (d = sqrt.(d))
    d = _additive_constant(d, add)

    # principal coordinates of -1/2 d² (double centred)
    x = d .^ 2
    for i in axes(x, 1)
        x[i, i] = 0.0
    end
    x = _dblcen(x)
    e = eigen(Symmetric(-x ./ 2))
    ord = sortperm(e.values; rev=true)
    vals = e.values[ord]
    want = abs.(vals) .> max(TOL, TOL * vals[1])
    eig = vals[want]
    vectors = e.vectors[:, ord[want]] .* sqrt.(abs.(eig))'
    pos = eig .> 0

    levels = sort(unique(g))
    idx = [findall(==(l), g) for l in levels]
    k = length(eig)
    cent = zeros(length(levels), k)
    for (a, ii) in enumerate(idx)
        if type === :centroid
            cent[a, :] = mean(vectors[ii, :]; dims=1)
        else
            cent[a, pos] = spatial_median(vectors[ii, pos])
            any(.!pos) && (cent[a, .!pos] = spatial_median(vectors[ii, .!pos]))
        end
    end
    gi = [searchsortedfirst(levels, x) for x in g]
    dpos = [sum(abs2, vectors[i, pos] .- cent[gi[i], pos]) for i in eachindex(g)]
    dneg = any(.!pos) ? [sum(abs2, vectors[i, .!pos] .- cent[gi[i], .!pos]) for i in eachindex(g)] :
           zeros(length(g))
    if any(dneg .> dpos)
        @warn "some squared distances are negative and changed to zero"
    end
    z = sqrt.(max.(dpos .- dneg, 0.0))
    if bias_adjust
        ng = length.(idx)
        z .*= sqrt.(ng[gi] ./ (ng[gi] .- 1))
    end
    gd = [mean(z[ii]) for ii in idx]
    return BetaDisper(z, g, levels, gd, cent, vectors, eig, type, bias_adjust, labs)
end

betadisper(X::AbstractMatrix, group::AbstractVector, metric::PreMetric; kw...) =
    betadisper(pairwise(metric, X; dims=1), group; kw...)
betadisper(X::AbstractMatrix, group::AbstractVector, metric::Type{<:PreMetric}; kw...) =
    betadisper(X, group, metric(); kw...)

"""
    dispersion(D, group; kw...)
    dispersion(X, group, metric; kw...)

The name of [`betadisper`](@ref) before version 0.2.
"""
dispersion(args...; kw...) = betadisper(args...; kw...)

function _dblcen(x)
    y = x .- mean(x; dims=1)
    return y .- mean(y; dims=2)
end

# vegan's addLingoes / addCailliez
function _additive_constant(d, add)
    (add === false || add === nothing) && return d
    add === true && (add = :lingoes)
    n = size(d, 1)
    off = .!Matrix{Bool}(I, n, n)
    if add === :lingoes
        G = _dblcen(-d .^ 2 ./ 2)
        ac = abs(min(0.0, minimum(eigvals(Symmetric(G)))))
        return ifelse.(off, sqrt.(d .^ 2 .+ 2ac), 0.0)
    elseif add === :cailliez
        G1 = _dblcen(-d .^ 2 ./ 2)
        G2 = _dblcen(-d ./ 2)
        M = [zeros(n, n) 2G1; -Matrix(1.0I, n, n) -4G2]
        ac = maximum(real.(eigvals(M)))
        return ifelse.(off, d .+ ac, 0.0)
    end
    throw(ArgumentError("add must be false, true, :lingoes or :cailliez, got $add"))
end
