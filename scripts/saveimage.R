save.image(paste0(".", snakemake@rule, ".RData"))


x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)
