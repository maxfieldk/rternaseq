log <- file(snakemake@log[[1]], open="wt")
sink(log)


library('DESeq2')
library('readr')
library('dplyr')
library("pheatmap")
library("ggplot2")
library('tibble')
library("RColorBrewer")

### inputs

counttypes = snakemake@params[["counttypes"]]
for (counttype in counttypes) {
    print(counttype)
    cts = read.delim(snakemake@input[[counttype]])
    rownames(cts) = cts$Geneid
    cts = select(cts, -Geneid)
    cnames = colnames(cts)

    coldata = read.csv(snakemake@params[["sampleinfo"]])

    contrasts = snakemake@params[["contrasts"]]

    levels = snakemake@params[["levels"]]

    outputdir = snakemake@params[["outputdir"]]
    ###

    condition = coldata$condition
    colnames(cts) <- coldata$sample
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ condition)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    #this sets prol as the reference level since its first in the vector
    dds$condition <- factor(dds$condition, levels = levels)

    ####
    dds <- DESeq(dds)
    ####
    resultsNames(dds) # lists the coefficients

    ####
    vsd <- vst(dds, blind=FALSE)
    sampleDists <- dist(t(assay(vsd)))

    pdf(paste(outputdir,counttype,"plots","pcaplot.pdf", sep = '/'), width=10, height=8)
    plotPCA(vsd, intgroup=c("condition"))
    dev.off()



    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    pdf(paste(outputdir,counttype,"plots","heatmapplot.pdf", sep = '/'), width=10, height=8)
    pheatmap(sampleDistMatrix,
            clustering_distance_rows=sampleDists, 
            clustering_distance_cols=sampleDists,
            col=colors)
    dev.off()
    ####


    for (contrast in contrasts) {
        res <- results(dds, name=contrast)
        # or to shrink log fold changes association with condition:
        resLFC <- lfcShrink(dds, coef=contrast, type="apeglm")
        # Save the plot as a PDF file
        pdf(paste(outputdir,counttype,contrast,"deplot.pdf", sep = '/'), width=10, height=8)
        plotMA(res, ylim=c(-7,7), main=paste0(contrast))
        dev.off()
        pdf(paste(outputdir,counttype,contrast,"deLFCplot.pdf", sep = '/'), width=10, height=8)
        plotMA(resLFC, ylim=c(-10,10), main=paste0(contrast))
        dev.off()
        

        resOrdered <- res[order(res$pvalue),]
        resSig <- subset(resOrdered, padj < 0.1)
        counttablesizenormed <- counts(dds, normalized=T)
        rld <- rlog(dds, blind=FALSE)

        write.csv(as.data.frame(resOrdered), file=paste(outputdir,counttype,contrast,"results.csv", sep='/'))
        write.csv(as.data.frame(resSig), file=paste(outputdir,counttype,contrast,"resultsSig.csv", sep='/'))
        write.csv(as.data.frame(counttablesizenormed), file=paste(outputdir,counttype,contrast,"counttablesizenormed.csv", sep='/'))
        write.csv(as.data.frame(assay(rld)), file=paste(outputdir,counttype,contrast,"rlogcounts.csv", sep='/'))

    }

}
