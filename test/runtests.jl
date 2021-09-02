using BetaDisp, Distances, StatsBase,BenchmarkTools
using Test

#@testset "BetaDisp.jl" begin
    x = rand(30,50) 
    y =rand(1:3,30)
    d = dispersion(x,y,Euclidean)
    permutest(d, 10000)
    dd = dispersion(x,y, metric = true)
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
