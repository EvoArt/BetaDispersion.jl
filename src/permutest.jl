"""
    DispersionTest

Result of [`permutest`](@ref) or [`dispersion_anova`](@ref): a one-way ANOVA
of the distances to group centres.

- `df`, `sumsq`, `meansq`: for groups and residuals
- `F`, `p`: the F statistic and its p-value (permutation or F distribution)
- `permutations`: number of permutations (0 for the parametric test)
- `pairwise`: `nothing`, or a named tuple of vectors `group1`, `group2`, `t`,
  `p` (two-sided permutation p-value), `p_adjusted` and `p_parametric`
  (pooled-variance t-test)
- `adjust`: the multiple-testing correction applied to `p_adjusted`
"""
struct DispersionTest
    df::Tuple{Int,Int}
    sumsq::Tuple{Float64,Float64}
    meansq::Tuple{Float64,Float64}
    F::Float64
    p::Float64
    permutations::Int
    pairwise::Union{Nothing,NamedTuple}
    adjust::Symbol
    group_distances::Vector{Float64}
    levels::Vector{String}
end

# one-way ANOVA of z on group index gi (k groups); returns (ssg, ssr, F)
function _anova(z, gi, k)
    n = length(z)
    s = zeros(k)
    c = zeros(Int, k)
    for i in 1:n
        s[gi[i]] += z[i]
        c[gi[i]] += 1
    end
    m = s ./ c
    mu = sum(z) / n
    ssg = sum(c[a] * (m[a] - mu)^2 for a in 1:k)
    ssr = sum((z[i] - m[gi[i]])^2 for i in 1:n)
    return ssg, ssr, (ssg / (k - 1)) / (ssr / (n - k))
end

function _pooled_t(x, y)
    m, n = length(x), length(y)
    pooled = sqrt(((m - 1) * var(x) + (n - 1) * var(y)) / (m + n - 2))
    return (mean(x) - mean(y)) / (pooled * sqrt(1 / m + 1 / n))
end

_group_index(bd) = [searchsortedfirst(bd.levels, g) for g in bd.group]

"""
    dispersion_anova(bd::BetaDisper)

Parametric one-way ANOVA of the distances to group centres, as
`anova(betadisper(...))` in vegan.
"""
function dispersion_anova(bd::BetaDisper)
    k = length(bd.levels)
    k >= 2 || throw(ArgumentError("need at least two groups"))
    n = length(bd.distances)
    ssg, ssr, F = _anova(bd.distances, _group_index(bd), k)
    df = (k - 1, n - k)
    p = fdistccdf(df[1], df[2], F)
    return DispersionTest(df, (ssg, ssr), (ssg / df[1], ssr / df[2]), F, p, 0, nothing, :none,
                          bd.group_distances, string.(bd.levels))
end

"""
    permutest(bd::BetaDisper; permutations=999, pairwise=false, adjust=:holm,
              rng=Random.default_rng())

Permutation test of differences in dispersion between groups. The F
statistic of a one-way ANOVA of the distances to group centres is compared
with its distribution when the ANOVA residuals are permuted (as vegan's
`permutest.betadisper`). P-values are `(1 + #{F_perm ≥ F}) / (1 + permutations)`.

With `pairwise=true` each pair of groups is also compared with a two-sided
permutation t-test (distances permuted between the two groups only), and
the p-values are corrected for multiple testing: `adjust` is `:holm`
(default), `:bh` (Benjamini-Hochberg), `:bonferroni` or `:none`.

`permutations` may also be a matrix with one permutation of `1:n` per row.
"""
function permutest(bd::BetaDisper; permutations=999, pairwise::Bool=false, adjust::Symbol=:holm,
                   rng::AbstractRNG=Random.default_rng())
    k = length(bd.levels)
    k >= 2 || throw(ArgumentError("need at least two groups"))
    z = bd.distances
    n = length(z)
    gi = _group_index(bd)
    ssg, ssr, F0 = _anova(z, gi, k)
    # residuals from the group means, permuted and refitted
    means = [mean(z[gi .== a]) for a in 1:k]
    res = z .- means[gi]
    perms = _perm_source(permutations, n)
    nperm = _nperm(perms)
    pr = similar(res)
    hits = 0
    EPS = sqrt(eps())
    for b in 1:nperm
        p = _perm(perms, b, rng, n)
        for i in 1:n
            pr[i] = res[p[i]]
        end
        _, _, Fb = _anova(pr, gi, k)
        hits += Fb >= F0 - EPS
    end
    pval = (hits + 1) / (nperm + 1)
    pw = pairwise ? _pairwise(bd, gi, k, nperm, rng, adjust) : nothing
    df = (k - 1, n - k)
    return DispersionTest(df, (ssg, ssr), (ssg / df[1], ssr / df[2]), F0, pval, nperm, pw, adjust,
                          bd.group_distances, string.(bd.levels))
end

function _pairwise(bd, gi, k, nperm, rng, adjust)
    z = bd.distances
    g1, g2, ts, ps, pp = String[], String[], Float64[], Float64[], Float64[]
    EPS = sqrt(eps())
    for a in 1:k-1, b in a+1:k
        x = z[gi .== a]
        y = z[gi .== b]
        t0 = _pooled_t(x, y)
        xy = vcat(x, y)
        m = length(x)
        hits = 0
        for _ in 1:nperm
            shuffle!(rng, xy)
            hits += abs(_pooled_t(view(xy, 1:m), view(xy, m+1:length(xy)))) >= abs(t0) - EPS
        end
        push!(g1, string(bd.levels[a]))
        push!(g2, string(bd.levels[b]))
        push!(ts, t0)
        push!(ps, (hits + 1) / (nperm + 1))
        push!(pp, 2 * tdistccdf(length(xy) - 2, abs(t0)))
    end
    return (group1=g1, group2=g2, t=ts, p=ps, p_adjusted=adjust_p(ps, adjust), p_parametric=pp)
end

"""
    adjust_p(p, method=:holm)

Adjust p-values for multiple testing: `:holm`, `:bh` (Benjamini-Hochberg),
`:bonferroni` or `:none`. The same as R's `p.adjust`.
"""
function adjust_p(p, method::Symbol=:holm)
    m = length(p)
    m == 0 && return Float64[]
    method === :none && return collect(float.(p))
    method === :bonferroni && return min.(1.0, m .* p)
    o = sortperm(p)
    q = similar(p, Float64)
    if method === :holm
        run = 0.0
        for (r, i) in enumerate(o)
            run = max(run, min(1.0, (m - r + 1) * p[i]))
            q[i] = run
        end
    elseif method === :bh
        run = 1.0
        for r in m:-1:1
            i = o[r]
            run = min(run, m / r * p[i])
            q[i] = min(run, 1.0)
        end
    else
        throw(ArgumentError("adjust must be :holm, :bh, :bonferroni or :none, got :$method"))
    end
    return q
end

# permutations: a count (random) or a matrix with one permutation per row
_perm_source(p::Integer, n) = p >= 0 ? Int(p) : throw(ArgumentError("permutations must be ≥ 0"))
function _perm_source(P::AbstractMatrix{<:Integer}, n)
    size(P, 2) == n || throw(DimensionMismatch("permutation matrix has $(size(P, 2)) columns, need $n"))
    return P
end
_nperm(p::Int) = p
_nperm(P::AbstractMatrix) = size(P, 1)
_perm(::Int, b, rng, n) = randperm(rng, n)
_perm(P::AbstractMatrix, b, rng, n) = view(P, b, :)
