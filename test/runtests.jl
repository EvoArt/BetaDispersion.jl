using BetaDisp, Distances, StatsBase,BenchmarkTools
using Test

#@testset "BetaDisp.jl" begin
    x = rand(100,50)
    y =rand(1:5,100)
    d = dispersion(x,y,Euclidean)
    permutest(d)
    
    bench = @benchmark dispersion($x,$y, BrayCurtis)
    @test mean(bench.times) .< 2
    bench = @benchmark permutest($d)
    @test mean(bench.times) .< 3
    
   distfuns =  [
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
         @test typeof(dispersion(x,rand(1:5,100), distfun)) <: Disp
    end

#end
