log <- file(snakemake@log[[1]], open="wt")
sink(log)

df = read.delim(snakemake@input[["telocal"]][1])
coldata = read.csv(snakemake@params[["sample_table"]])

bounddf = data.frame(Geneid = df[,1])
for (sample in snakemake@input[["telocal"]]) {
    bounddf = cbind(bounddf, read.delim(sample)[,2])
}
colnames(bounddf) = c("Geneid", coldata$sample_name)
write.table(bounddf, file = snakemake@output[["aggcounts"]], row.names = FALSE, sep = '\t')
