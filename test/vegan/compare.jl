# Compare BetaDispersion.jl with vegan::betadisper on reference.R's output.
# Usage: julia --project=. test/vegan/compare.jl OUTDIR
using BetaDispersion, Distances, DelimitedFiles, Test
dir = ARGS[1]
rd(f) = readdlm(joinpath(dir, f), ','; header=true)[1]
g = String.(rd("group.csv")[:, 1])
Y = Float64.(rd("Y.csv"))
P = Int.(rd("perm.csv"))
D = pairwise(BrayCurtis(), Y; dims=1)
@testset "vegan betadisper" begin
    for (type, tol) in ((:centroid, 1e-10), (:median, 1e-5))
        @testset "$type" begin
            bd = betadisper(D, g; type=type)
            s = rd("summary_$type.csv")
            # vegan's median comes from optim(BFGS) with reltol 1e-8, so single
            # distances differ by up to ~1e-4; group means agree to < 1e-6
            @test bd.distances ≈ Float64.(rd("dist_$type.csv")[:, 1]) rtol = (type === :median ? 1e-3 : tol)
            @test bd.levels == String.(s[:, 1])
            @test bd.group_distances ≈ Float64.(s[:, 2]) rtol = tol
            t = permutest(bd; permutations=P)
            @test t.F ≈ s[1, 3] rtol = 100tol
            type === :centroid && @test t.p == s[1, 4]
        end
    end
    @test betadisper(D, g; add=:lingoes).distances ≈ Float64.(rd("dist_lingoes.csv")[:, 1]) rtol = 1e-3
    @test betadisper(D, g; add=:cailliez).distances ≈ Float64.(rd("dist_cailliez.csv")[:, 1]) rtol = 1e-3
end
