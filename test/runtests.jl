using BetaDispersion, Distances, StatsBase,BenchmarkTools, RCall
using Test
R"library(vegan)"

   x = rand(1000,100) 
   y =rand(1:20,1000)
   j = @benchmark dispersion($x,$y,BrayCurtis) 
   r = @benchmark R"betadisper(vegdist($x),as.factor($y))" 
   @test mean(j.times) < mean(r.times)
   
   Dj = pairwise(BrayCurtis(),x,dims = 1) 
   Dr =R"vegdist($x)"   
   dispj = dispersion(Dj,y)
   dispr = R"betadisper($Dr,as.factor($y))"
   j = @benchmark permutest(dispj,9999) 
   r = @benchmark R"permutest($dispr,pairwise = TRUE, permutations = 9999)" 
   @test mean(j.times) < mean(r.times)

