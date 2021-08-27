using BetaDisp, Distances, StatsBase,BenchmarkTools
using Test

#@testset "BetaDisp.jl" begin
    x = rand(30,50) .*y
    y =rand(1:5,30)
    d = dispersion(x,y,Euclidean)
    permutest(d)
    
    bench = @benchmark dispersion($x,$y, BrayCurtis)
    @test mean(bench.times) .< 2
    bench = @benchmark permutest($d)
    @test mean(bench.times) .< 3
    
   distfuns =  [
    Euclidean,
    Cityblock,
    TotalVariation,
    Chebyshev,
    Minkowski,
    Jaccard,
    BrayCurtis]
    for distfun in distfuns
         @test typeof(dispersion(x,rand(1:5,100), distfun)) <: Disp
    end

#end
