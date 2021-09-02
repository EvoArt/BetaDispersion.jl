using BetaDispersion, Distances, StatsBase,BenchmarkTools
using Test

@testset "BetaDispersion.jl" begin
    x = rand(1000,50) 
    y =rand(1:10,1000)
    d = dispersion(x,y,Euclidean)
    dispersion(x,y, metric = true)
    bench = @benchmark dispersion($x,$y, BrayCurtis)
    @test mean(bench.times) < 2*10^9
    bench = @benchmark permutest($d)
    @test mean(bench.times) < 4*10^9
    
end
