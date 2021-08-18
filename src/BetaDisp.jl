module BetaDisp

using LinearAlgebra,  Random, DirectionalStatistics, NamedArrays
include("ANOVA.jl")
include("dispersion.jl")
include("permutest.jl")

export dispersion, permutest
end
