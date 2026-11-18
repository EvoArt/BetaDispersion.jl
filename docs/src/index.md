```@meta
CurrentModule = BetaDispersion
```

# BetaDispersion.jl

`BetaDispersion` tests whether groups differ in their spread in a dissimilarity
space. It follows `vegan::betadisper` and its permutation test. A significant
result can help distinguish differences in dispersion from differences in group
location when interpreting PERMANOVA.

## Quick start

```julia
using BetaDispersion

# D is a square dissimilarity matrix; group has one label per row of D.
bd = betadisper(D, group)
permutest(bd; permutations=999, pairwise=true)
dispersion_anova(bd)
```

Pass observations and a `Distances.jl` metric instead with
`betadisper(X, group, BrayCurtis())`, where rows of `X` are observations. The
default group centre is the spatial median; use `type=:centroid` for centroids.

## API

```@docs
BetaDispersion
BetaDisper
DispersionTest
betadisper
permutest
dispersion_anova
adjust_p
dispersion
```

```@index
```
