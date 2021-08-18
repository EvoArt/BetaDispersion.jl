using BetaDisp, Distances, StatsFuns
using Test

@testset "BetaDisp.jl" begin
    x = rand(1000,50)
    d = dispersion(x,rand(1:5,1000))
    bench = @benchmark dispersion(x,rand(1:5,1000), BrayCurtis)
    @test mean(bench.times) .< 0.1
    bench = @benchmark permutest(d)
    @test mean(bench.times) .< 1
end
