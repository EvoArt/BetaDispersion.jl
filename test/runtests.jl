using BetaDisp, Distances, StatsBase,BenchmarkTools
using Test

#@testset "BetaDisp.jl" begin
    x = rand(1000,50)
    d = dispersion(x,rand(1:5,1000),BrayCurtis)
    bench = @benchmark dispersion(x,rand(1:5,1000), BrayCurtis)
    @test mean(bench.times) .< 0.1
    bench = @benchmark permutest(d)
    @test mean(bench.times) .< 1
    
   distfuns =  [Euclidean,
    SqEuclidean,
    PeriodicEuclidean,
    Cityblock,
    TotalVariation,
    Chebyshev,
    Minkowski,
    Jaccard,
    BrayCurtis,
    RogersTanimoto,
    Hamming,
    CosineDist,
    CorrDist,
    ChiSqDist,
    KLDivergence,
    GenKLDivergence,
    JSDivergence,
    RenyiDivergence,
    SpanNormDist,
    WeightedEuclidean,
    WeightedSqEuclidean,
    WeightedCityblock,
    WeightedMinkowski,
    WeightedHamming,
    SqMahalanobis,
    Mahalanobis,
    BhattacharyyaDist,
    HellingerDist,
    Haversine,
    SphericalAngle,
    MeanAbsDeviation,
    MeanSqDeviation,
    RMSDeviation,
    NormRMSDeviation,
    Bregman]
    for distfun in distfuns
        @test dispersion(x,rand(1:5,1000), distfun)
    end

#end
