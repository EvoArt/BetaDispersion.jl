# BetaDispersion

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://EvoArt.github.io/BetaDispersion.jl/dev/)
[![CI](https://github.com/EvoArt/BetaDispersion.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/EvoArt/BetaDispersion.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/EvoArt/BetaDispersion.jl/graph/badge.svg)](https://codecov.io/gh/EvoArt/BetaDispersion.jl)

Multivariate homogeneity of group dispersions (Anderson 2006), as
`vegan::betadisper`, with the permutation test of `vegan::permutest`. Often
used in ecology to compare beta diversity between groups of sites.

```julia
using Pkg
Pkg.add("BetaDispersion")
```

```julia
using BetaDispersion, Distances

bd = betadisper(X, group, BrayCurtis())   # rows of X are observations; or betadisper(D, group)
permutest(bd; pairwise=true)              # global F test and pairwise t tests
dispersion_anova(bd)                      # parametric ANOVA
```

Distances to group centres use spatial medians by default (`type=:centroid`
for centroids). Pairwise p-values are adjusted with `adjust_p` (Holm by
default).

Anderson MJ (2006). Distance-based tests for homogeneity of multivariate
dispersions. *Biometrics* 62, 245–253.
[doi:10.1111/j.1541-0420.2005.00440.x](https://doi.org/10.1111/j.1541-0420.2005.00440.x)
