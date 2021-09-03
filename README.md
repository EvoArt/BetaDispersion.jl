# BetaDispersion


This is a small package aimed at providing equivalent functionality to `betadisper` in the R package `vegan`, namely multivariate test for homogeneity of variance (dispersion). This is often used in ecology to compare beta diversity between metacommunities. This analysis comes from work by [Marti Anderson (2006)](https://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2005.00440.x). However, a number of tweaks have been made to the preferred way of performing the analysis since then. We aim to keep this package aligned with decisions made by [vegan devs](https://github.com/vegandevs/vegan/blob/master/R/betadisper.R). However, the current implementation only offers permutation tests for P-values, and only spatial medians (not centroids) are used. Use of permutation tests and spatial medians are considered more robust than alternatives, but please open an issue/PR if you feel that these options are not appropriate for your use case. 

The basic steps are to calculate spatial medians for groups of data points, calculate F-statistics (including those for each pairwise combination of groups), then permute the residuals (distances from medians) and re-calculate F for each permutation. The P-value is the proportion of permuted F values that 
are greater than or or equal to the original F value.

<img src="https://github.com/EvoArt/BetaDispersion.jl/blob/master/docs/permutation.gif"/>

Two functions are exported, `dispersion` takes either a data array (where each row is an observation), a vector containing group identities and a distance function type from `Distances.jl`, or a distance matrix and the grouping vector. Alternatively, if working with Euclidean distances, pass in the original data (where each row is an observation) instead of a distance matrix and set the key work argument `metric = true`. This function returns a `Disp` struct containing:
*    `F` = Global F-Statistic 
*    `pairwise_F` = pairwise F-statistics
*    `medians` = spatial (aka geometric) median of each group in the transformed coordinates
*    `residuals` = vector containing each observations distance from its group median
*    `means` = vector containing the means distance to median for each group
*    `group` = the original grouping vector
*    `levels` = unique items from 'group'

`permutest` takes the named tuple returned by `dispersion` and returns a named tuple containing:
 *   `P` = Global P-value
 *   `pairwise_P` = pairwise P-values
 *   `F` = Global F-Statistic 
 *   `pairwise_F` = pairwise F-statistics

Additionally, If [Turing.jl](https://turing.ml/stable/) is loaded before BetaDispersion you can call `bayesdisp` on the `Disp` returned by `dispersion` to perform Bayesian inference to estimate mean and standard deviation of the residuals of each group. By default 4 chains are run simultaneously using the No U-turns (NUTS) sampler for 1000 iterations each. The function returns the Markov chains along with diagnostics, more details can be found [here](https://turinglang.github.io/MCMCChains.jl/dev/). Check the [source code](https://github.com/EvoArt/BetaDispersion.jl/blob/master/src/Bayes.jl) to see the function arguments. If fine control over model implementation is required, the user is urged to construct their model directly in Turing.jl or the users preferred PPL.

## Example usage

```julia
using Turing, StatsPlots
using Distances, BetaDispersion
x = rand(30,5)
g = rand(1:2,30)
d = dispersion(x,g,BrayCurtis)
p = permutest(d)
chns = bayesdisp(d)
plot(chns)
```
<img src="https://github.com/EvoArt/BetaDispersion.jl/blob/master/docs/example.png">


