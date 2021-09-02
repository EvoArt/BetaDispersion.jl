using BetaDisp, Distances, StatsBase,BenchmarkTools
using Test

@testset "BetaDisp.jl" begin
    x = rand(500,50) 
    y =rand(1:10,500)
    d = dispersion(x,y,Euclidean)
    dispersion(x,y, metric = true)
    bench = @benchmark dispersion($x,$y, BrayCurtis)
    @test mean(bench.times) < 1*10^8
    bench = @benchmark permutest($d)
    @test mean(bench.times) < 1*10^9
    
end
