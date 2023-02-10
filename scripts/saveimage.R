save.image()
x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)
