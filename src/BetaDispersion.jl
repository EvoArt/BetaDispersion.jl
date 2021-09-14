module BetaDispersion

using LinearAlgebra,  Random, DirectionalStatistics, NamedArrays, Distances, Requires
import StatsBase: mean
include("ANOVA.jl")
include("dispersion.jl")
include("permutest.jl")
function __init__()
    @require Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0" include("Bayes.jl")
end

export dispersion, permutest
end
