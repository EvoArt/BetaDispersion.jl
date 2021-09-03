using BetaDispersion, Distances, StatsBase,BenchmarkTools, RCall
using Test
vegan = "vegan"
R"install.packages($vegan)"
R"library(vegan)"
# Rcall uses a single thread, so make sure Julia is using only ine thread as well
@testset "BetaDispersion.jl" begin
    #test against vegan on small data set
    x = rand(30,10) 
    y =rand(1:3,30)
    Dj = pairwise(BrayCurtis(),x,dims = 1)
    Dr =R"vegdist($x)"
    j = @benchmark dispersion(Dj,y)
    r = @benchmark R"betadisper($Dr,as.factor($y))"
    @test mean(j.times) < mean(r.times)
    dispj = dispersion(Dj,y)
    dispr = R"betadisper($Dr,as.factor($y))"
    j = @benchmark permutest(dispj,999)
    r = @benchmark R"permutest($dispr,pairwise = TRUE, permutations = 999)"
    @test mean(j.times) < mean(r.times)

    #try with more variables/species
    x = rand(30,100) 
    y =rand(1:3,30)
    Dj = pairwise(BrayCurtis(),x,dims = 1)
    Dr =R"vegdist($x)"
    j = @benchmark dispersion(Dj,y)
    r = @benchmark R"betadisper($Dr,as.factor($y))"
    @test mean(j.times) < mean(r.times)
    dispj = dispersion(Dj,y)
    dispr = R"betadisper($Dr,as.factor($y))"
    j = @benchmark permutest(dispj,999) 
    r = @benchmark R"permutest($dispr,pairwise = TRUE, permutations = 999)" 
    @test mean(j.times) < mean(r.times)

   #now with a larger data set and more groups  
   x = rand(1000,100) 
   y =rand(1:20,1000)
   Dj = pairwise(BrayCurtis(),x,dims = 1) #  Distances.jl returns a distance matrix insantaneously
   Dr =R"vegdist($x)"   # vegan takes a very long time.
   j = @benchmark dispersion(Dj,y) #318.505 ms
   r = @benchmark R"betadisper($Dr,as.factor($y))" #4.477 s
   @test mean(j.times) < mean(r.times)
   dispj = dispersion(Dj,y)
   dispr = R"betadisper($Dr,as.factor($y))"
   j = @benchmark permutest(dispj,999) #279.505 ms on my machine (single thread)
   r = @benchmark R"permutest($dispr,pairwise = TRUE, permutations = 999)" #11.515 s
   @test mean(j.times) < mean(r.times)


  
end
