"""
    BetaDispersion

Multivariate homogeneity of group dispersions (Anderson 2006), as
`vegan::betadisper`, with the permutation test of `vegan::permutest`.

```julia
bd = betadisper(D, group)             # D: a dissimilarity matrix
permutest(bd; pairwise=true)
```
"""
module BetaDispersion

using LinearAlgebra, Random, Statistics
using Distances: pairwise, PreMetric
using Printf: @sprintf
using StatsFuns: fdistccdf, tdistccdf

export betadisper, permutest, dispersion_anova, adjust_p, BetaDisper, DispersionTest
export dispersion   # name used before 0.2

include("median.jl")
include("betadisper.jl")
include("permutest.jl")
include("show.jl")

end
