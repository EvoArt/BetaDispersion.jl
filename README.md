# BetaDisp

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://EvoArt.github.io/BetaDisp.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://EvoArt.github.io/BetaDisp.jl/dev)
[![Build Status](https://github.com/EvoArt/BetaDisp.jl/workflows/CI/badge.svg)](https://github.com/EvoArt/BetaDisp.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/EvoArt/BetaDisp.jl?svg=true)](https://ci.appveyor.com/project/EvoArt/BetaDisp-jl)
[![Coverage](https://codecov.io/gh/EvoArt/BetaDisp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/EvoArt/BetaDisp.jl)

This is a small package aiming at providing equivalent functionality to `betadisper` in the R package `vegan`, namely multivariate test for homogeneity of variance (dispersion). This is often used in ecology to compare beta diversity between metacommunities. This analysis comes from work by [Maritin Anderson (2006)](https://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2005.00440.x). However, a number of tweaks have been made to the preferred way of performing the analyis since then. We aim to keep this package aligned with decisions made by [vegan devs](https://github.com/vegandevs/vegan/blob/master/R/betadisper.R). However, the current implementation only offers permutation tests for P-values, and only spatial medians (not centroids) are used. 

Two functions are exported, `dispersion` takes either an array, a vector containing group identities and a distance function type from `Distances.jl`, or a distance matrix a and the grouping vector. This returns a named tuple containing:
*    F = Global F-Statistic 
*    pairwise_F = pairwise F-statistics
*    medians = spatial (aka geometric) median of each group in the transformed coordinates
*    residuals = vector containing each observations distance from its group median
*    means = vector containing the means distance to median for each group
*    group = the original grouping vector
*    levels = unique items from 'group'

`permutest` takes the named tuple returned by `dispersion` and returns a named tuple containing:
 *   P = Global P-value
 *   pairwise_P = pairwise P-value
 *   F = Global F-Statistic 
 *   pairwise_F = pairwise F-statistics

## Example usage

```julia
using Distance, BetaDisp
x = rand(9,5)
g = rand(1:2,9)
d = dispersion(x,g,BrayCurtis)
p = permutest(d)
```
