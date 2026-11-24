using BetaDispersion, Test, Distances, Random, Statistics, LinearAlgebra, StableRNGs
using BetaDispersion: spatial_median

sumdist(P, y) = sum(norm(P[i, :] .- y) for i in axes(P, 1))

@testset "BetaDispersion" begin
    @testset "spatial median" begin
        rng = StableRNG(1)
        for k in (1, 2, 5)
            P = randn(rng, 15, k)
            y = spatial_median(P)
            f = sumdist(P, y)
            for _ in 1:50
                @test f <= sumdist(P, y .+ 1e-4 .* randn(rng, k)) + 1e-12
            end
        end
        # median on a data point (Vardi-Zhang case): a cluster plus outliers
        P = [0.0 0.0; 0.0 0.0; 0.0 0.0; 1.0 0.0; 0.0 1.0; -1.0 0.0]
        @test norm(spatial_median(P)) < 1e-8
        @test spatial_median(reshape([2.0, 4.0, 9.0], 3, 1)) ≈ [4.0]
    end

    @testset "Euclidean oracle" begin
        rng = StableRNG(2)
        X = randn(rng, 40, 4)
        g = rand(rng, ["c", "a", "b"], 40)          # unsorted labels: the 0.1 bug
        D = pairwise(Euclidean(), X; dims=1)
        for type in (:centroid, :median)
            bd = betadisper(D, g; type=type)
            @test bd.levels == ["a", "b", "c"]
            for (a, l) in enumerate(bd.levels)
                ii = findall(==(l), g)
                c = type === :centroid ? vec(mean(X[ii, :]; dims=1)) : spatial_median(X[ii, :])
                @test bd.distances[ii] ≈ [norm(X[i, :] .- c) for i in ii] rtol = 1e-7
                @test bd.group_distances[a] ≈ mean(bd.distances[ii])
            end
        end
        # the result for an observation does not depend on the order of the data
        p = randperm(rng, 40)
        b1 = betadisper(D, g)
        b2 = betadisper(D[p, p], g[p])
        @test b2.distances ≈ b1.distances[p] rtol = 1e-7
        # data + metric form
        @test betadisper(X, g, Euclidean()).distances ≈ b1.distances
        @test dispersion(X, g, Euclidean).distances ≈ b1.distances
    end

    @testset "non-Euclidean dissimilarities, options" begin
        rng = StableRNG(3)
        Y = rand(rng, 0:20, 30, 12)
        g = repeat(1:3, 10)
        D = pairwise(BrayCurtis(), Y; dims=1)
        bd = betadisper(D, g)
        @test any(<(0), bd.eig)                    # Bray-Curtis is not Euclidean
        @test all(>=(0), bd.distances)
        ba = betadisper(D, g; bias_adjust=true)
        @test ba.distances ≈ bd.distances .* sqrt(10 / 9)
        for add in (:lingoes, :cailliez)
            b = betadisper(D, g; add=add)
            @test all(>=(-1e-8), b.eig)            # corrected to Euclidean
        end
        @test all(>(0), betadisper(D, g; sqrt_dist=true).distances)
        gm = Vector{Union{Missing,Int}}(g); gm[1] = missing
        @test length(betadisper(D, gm).distances) == 29
    end

    @testset "tests" begin
        rng = StableRNG(4)
        # group 2 is much more dispersed
        X = vcat(0.1 .* randn(rng, 20, 3), 2.0 .* randn(rng, 20, 3), 0.1 .* randn(rng, 20, 3))
        g = repeat(["a", "b", "c"]; inner=20)
        bd = betadisper(pairwise(Euclidean(), X; dims=1), g)
        t = permutest(bd; pairwise=true, permutations=499, rng=StableRNG(5))
        @test t.p == 1 / 500
        @test t.permutations == 499
        a = dispersion_anova(bd)
        @test a.F ≈ t.F
        @test a.p < 1e-6
        pw = t.pairwise
        @test pw.group1 == ["a", "a", "b"] && pw.group2 == ["b", "c", "c"]
        @test pw.p[1] < 0.01 && pw.p[3] < 0.01 && pw.p[2] > 0.01
        @test all(pw.p_adjusted .>= pw.p)
        # no difference: p-values are roughly uniform
        Z = randn(rng, 30, 3)
        bz = betadisper(pairwise(Euclidean(), Z; dims=1), repeat(1:2, 15))
        @test permutest(bz; permutations=199, rng=StableRNG(6)).p > 0.005
        # a permutation matrix fixes the result
        P = reduce(vcat, (randperm(StableRNG(i), 60)' for i in 1:99))
        @test permutest(bd; permutations=P).p == permutest(bd; permutations=P).p
        @test occursin("Pr(>F)", sprint(show, MIME("text/plain"), t))
        @test occursin("spatial median", sprint(show, MIME("text/plain"), bd))
    end

    @testset "p-value adjustment (R p.adjust)" begin
        p = [0.01, 0.04, 0.03, 0.2]
        @test adjust_p(p, :holm) ≈ [0.04, 0.09, 0.09, 0.2]
        @test adjust_p(p, :bh) ≈ [0.04, 0.16 / 3, 0.16 / 3, 0.2]
        @test adjust_p(p, :bonferroni) ≈ [0.04, 0.16, 0.12, 0.8]
        @test adjust_p(p, :none) == p
        @test_throws ArgumentError adjust_p(p, :fdr)
    end

    @testset "edge cases and displays" begin
        # Weiszfeld iterate starts on a data point that is not the median
        P = [0.0 0.0; 0.0 1.0; 0.0 -1.0; 3.0 0.0; 3.0 0.5]
        y = spatial_median(P)
        @test norm(y) > 0.1
        rng = StableRNG(7)
        for _ in 1:50
            @test sumdist(P, y) <= sumdist(P, y .+ 1e-4 .* randn(rng, 2)) + 1e-12
        end
        # strongly non-Euclidean dissimilarities: negative squared distances
        D = [0.0 0.381 0.906 0.496 0.316 0.045; 0.381 0.0 0.001 0.019 1.254 0.077;
             0.906 0.001 0.0 0.816 0.25 0.082; 0.496 0.019 0.816 0.0 1.136 0.046;
             0.316 1.254 0.25 1.136 0.0 0.963; 0.045 0.077 0.082 0.046 0.963 0.0]
        bd = @test_logs (:warn, r"negative") match_mode = :any betadisper(D, [1, 1, 1, 2, 2, 2])
        @test all(>=(0), bd.distances)
        @test_throws ArgumentError betadisper(D, [1, 1, 1, 2, 2, 2]; add=:other)
        @test_throws ArgumentError betadisper(D, [1, 1, 1, 2, 2, 2]; type=:mean)
        # compact displays
        @test sprint(show, bd) == "BetaDisper(6 observations, 2 groups)"
        t = permutest(bd; permutations=9, rng=StableRNG(8))
        @test startswith(sprint(show, t), "DispersionTest(F = ")
        @test occursin("centroid", sprint(show, MIME("text/plain"), betadisper(D, [1, 1, 1, 2, 2, 2]; type=:centroid)))
    end
end
