# Reference results from vegan::betadisper for test/vegan/compare.jl.
# Usage: Rscript reference.R OUTDIR
suppressMessages(library(vegan))
out <- commandArgs(TRUE)[1]
set.seed(2)
n <- 60
g <- sample(c("b", "a", "c"), n, TRUE)
Y <- matrix(rpois(n * 8, 6), n, 8)
Y[g == "c", ] <- Y[g == "c", ] + matrix(rpois(sum(g == "c") * 8, 6), ncol = 8)
D <- vegdist(Y, "bray")
perm <- shuffleSet(n, 999)
w <- function(x, f) write.csv(x, file.path(out, f), row.names = FALSE)
w(data.frame(g = g), "group.csv"); w(Y, "Y.csv"); w(unclass(perm), "perm.csv")
for (type in c("median", "centroid")) {
  bd <- betadisper(D, g, type = type)
  pt <- permutest(bd, permutations = perm)
  w(data.frame(d = bd$distances), paste0("dist_", type, ".csv"))
  w(data.frame(level = names(bd$group.distances), gd = bd$group.distances,
               F = pt$tab$F[1], P = pt$tab[["Pr(>F)"]][1]), paste0("summary_", type, ".csv"))
}
bl <- betadisper(D, g, add = "lingoes"); w(data.frame(d = bl$distances), "dist_lingoes.csv")
bc <- betadisper(D, g, add = "cailliez"); w(data.frame(d = bc$distances), "dist_cailliez.csv")
cat("vegan", as.character(packageVersion("vegan")), "\n")
