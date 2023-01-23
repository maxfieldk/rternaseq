##multi
df = read.delim(snakemake@input[["telocalMULTI"]][1])
coldata = read.csv(snakemake@params[["sample_table"]])

bounddf = data.frame(Geneid = df[,1])
for (sample in snakemake@input[["telocal"]]) {
    bounddf = cbind(bounddf, read.delim(sample)[,2])
}
colnames(bounddf) = c("Geneid", coldata$sample_name)
write.table(bounddf, file = snakemake@output[["aggcountsMULTI"]], row.names = FALSE, sep = '\t')

##uniq
df = read.delim(snakemake@input[["telocalUNIQUE"]][1])
coldata = read.csv(snakemake@params[["sample_table"]])

bounddf = data.frame(Geneid = df[,1])
for (sample in snakemake@input[["telocal"]]) {
    bounddf = cbind(bounddf, read.delim(sample)[,2])
}
colnames(bounddf) = c("Geneid", coldata$sample)
write.table(bounddf, file = snakemake@output[["aggcountsUNIQUE"]], row.names = FALSE, sep = '\t')