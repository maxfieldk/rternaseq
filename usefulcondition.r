
rule featureCounts:
    input:
        sortedSTARbams = expand("../senescence/outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = samples) + expand("../arna/outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = samples)
    output:
        countsmessy = "outs/agg/refseq.counts_messy.txt",
        counts = "outs/agg/refseq.counts.txt",
        metafeaturecounts = "outs/agg/refseq.metafeature.counts.txt"
    params: 
        gtf = config['refseq']
    log: "logs/agg/featureCounts.log"
    conda:
        "envs/deeptools.yml"
    threads: 2
    shell: 
        """
featureCounts -p -T {threads} -t exon -a {params.gtf} -o {output.countsmessy} {input.sortedSTARbams} 2> {log}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
featureCounts -p -T {threads} -B -O -a {params.gtf} -o {output.metafeaturecounts} {input.sortedSTARbams} 2>> {log}
        """

rule mergeTElocal:
    input:
        telocalMULTIsource1 = expand("../senescence/outs/{sample}/TElocal/{sample}_multi.cntTable", sample = samples),
        telocalMULTIsource2 = expand("../arna/outs/{sample}/TElocal/{sample}_multi.cntTable", sample = samples),

        telocalUNIQUEsource1 = expand("../senescence/outs/{sample}/TElocal/{sample}_uniq.cntTable", sample = samples)
        telocalUNIQUEsource2 = expand("../arna/outs/{sample}/TElocal/{sample}_uniq.cntTable", sample = samples)
    params:
        sampleinfo = config["sample_table"]
    conda: "envs/renv.yml"
    log: "logs/agg/mergeTElocal.log"
    output:
        aggcountsMULTI = "outs/agg/TElocalCounts_MULTI.txt",
        aggcountsUNIQUE = "outs/agg/TElocalCounts_UNIQUE.txt"
    script:
        "scripts/mergeTEdf.R"    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############# in deseq sript
    condition = coldata$condition
    colnames(cts) <- coldata$sample_name
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ batch + condition)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    #this sets prol as the reference level since its first in the vector
    dds$condition <- factor(dds$condition, levels = levels)

    ####
    dds <- DESeq(dds)
    ####
    resultsNames(dds) # lists the coefficients

    ####

    vsd <- vst(dds)
    
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    
    pdf(paste(outputdir,counttype,"plots","pcaplot.pdf", sep = '/'), width=10, height=8)
    print(plotPCA(vsd, intgroup=c("condition", "batch")) + theme_cowplot())
    dev.off()

    assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
    plotPCA(vsd, "batch")
