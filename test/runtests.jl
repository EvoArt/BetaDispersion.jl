using BetaDispersion, Distances, StatsBase,BenchmarkTools, RCall
using Test
R"library(vegan)"

   x = rand(1000,100) 
   y =rand(1:20,1000)
   Dj = pairwise(BrayCurtis(),x,dims = 1) 
   Dr =R"vegdist($x)"   
   j = @benchmark dispersion(Dj,y) 
   r = @benchmark R"betadisper($Dr,as.factor($y))" 
   @test mean(j.times) < mean(r.times)
   dispj = dispersion(Dj,y)
   dispr = R"betadisper($Dr,as.factor($y))"
   j = @benchmark permutest(dispj,9999) 
   r = @benchmark R"permutest($dispr,pairwise = TRUE, permutations = 9999)" 
   @test mean(j.times) < mean(r.times)

