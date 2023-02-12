library(bedr)
index <- get.example.regions()
index[["a"]][1]

path <- "/users/mkelsey/data/ref/genomes/hs1/annotations2/repeatMaskerAnnotated.bed"

a <- "chr1:17011-17323"
a.int2 <- bedr(input = list(a = a, b = path), method = "intersect", params = "-loj -sorted")
