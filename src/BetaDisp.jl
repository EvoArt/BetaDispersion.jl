module BetaDisp

using LinearAlgebra,  Random, DirectionalStatistics, NamedArrays
include("src\\ANOVA.jl")
include("src\\dispersion.jl")
include("src\\permutest.jl")

export dispersion, permutest
end
