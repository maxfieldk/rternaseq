from snakemake.shell import shell
shell.executable('/usr/bin/bash')
import os
from pathlib import Path
configfile: "config.yaml"
pepfile: "conf/project_config.yaml"




container: "docker://maxfieldkelsey/snakesenescence:latest"

samples = pep.sample_table.sample_name

rule all:
    input:
        results = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"])
        


# tips:
# build needed directories in python! Will save you much frustration :)
# watch out for multiline shell commands - they fail for weird reasons sometimes.

########################################################### Make folders
if True:
    if not os.path.exists('envs'):
        os.makedirs('envs')
    if not os.path.exists('temp'):
        os.makedirs('temp')

    if not os.path.exists('results'):
        os.makedirs('results')

    if not os.path.exists('report'):
        os.makedirs('report')

    if not os.path.exists('results/agg'):
        os.makedirs('results/agg')

    if not os.path.exists('results/agg/deseq2'):
        os.makedirs('results/agg/deseq2')

    path = "results/agg/repeatanalysis"
    if not os.path.exists(path):
        os.makedirs(path)

        for telocaltype in config["telocaltypes"]:
            path = "results/agg/repeatanalysis/%s"%(telocaltype)
            if not os.path.exists(path):
                os.makedirs(path)

        for telocaltype in config["telocaltypes"]:
            for contrast in config["contrasts"]:
                path = 'results/agg/repeatanalysis/%s/%s'%(telocaltype,contrast)
                if not os.path.exists(path):
                    os.makedirs(path)

    for counttype in config["counttypes"]:
        if not os.path.exists('results/agg/deseq2/%s'%(counttype)):
            os.makedirs('results/agg/deseq2/%s'%(counttype))
        if not os.path.exists('results/agg/deseq2/%s/plots'%(counttype)):
            os.makedirs('results/agg/deseq2/%s/plots'%(counttype))
        
        for contrast in config["contrasts"]:
            if not os.path.exists('results/agg/deseq2/%s/%s'%(counttype, contrast)):
                os.makedirs('results/agg/deseq2/%s/%s'%(counttype, contrast))
                
    if not os.path.exists('results/agg/clusterprofiler'):
        os.makedirs('results/agg/clusterprofiler')

    for contrast in config["contrasts"]:
        if not os.path.exists('results/agg/clusterprofiler/%s'%(contrast)):
            os.makedirs('results/agg/clusterprofiler/%s'%(contrast))

    for contrast in config["contrasts"]:
        path = 'results/agg/clusterprofiler/%s/gsea'%(contrast)
        if not os.path.exists(path):
            os.makedirs(path)

    for contrast in config["contrasts"]:
        path = 'results/agg/clusterprofiler/%s/hypgeo'%(contrast)
        if not os.path.exists(path):
            os.makedirs(path)

        for contrast in config["contrasts"]:
            path = 'results/agg/clusterprofiler/%s/KEGG'%(contrast)
            if not os.path.exists(path):
                os.makedirs(path)

    if not os.path.exists('results/agg/plots'):
        os.makedirs('results/agg/plots')

    if not os.path.exists('rawdata'):
        os.makedirs('rawdata')

    if not os.path.exists('outs'):
        os.makedirs('outs')

    if not os.path.exists('outs/agg'):
        os.makedirs('outs/agg')

    if not os.path.exists('logs'):
        os.makedirs('logs')

    if not os.path.exists('logs/agg'):
        os.makedirs('logs/agg')

    if not os.path.exists('fastqcOut'):
        os.makedirs('fastqcOut')

    for e in samples:
        respath = 'results/%s' % (e)
        if not os.path.exists(respath):
            os.makedirs(respath)

        resplotpath = 'results/%s/plots' % (e)
        if not os.path.exists(resplotpath):
            os.makedirs(resplotpath)

        outpath = 'outs/%s' % (e)
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        
        outpath = 'outs/%s/TElocal' % (e)
        if not os.path.exists(outpath):
            os.makedirs(outpath)

        outSTARpath = 'outs/%s/star_output' % (e)
        if not os.path.exists(outSTARpath):
            os.makedirs(outSTARpath)
            Path(outSTARpath).chmod(0o0777)
        logspath = 'logs/%s' % (e)
        if not os.path.exists(logspath):
            os.makedirs(logspath)


###########################################################

rule prefetch:
    output:
        sra = temp("rawdata/{sample}/{sample}.sra")
    threads: 2
    conda:
        "envs/deeptools.yml"
    shell: "prefetch {wildcards.sample} --output-directory rawdata"

rule fastqdump:
    input: "rawdata/{sample}/{sample}.sra"
    threads: 2
    params:
        outdir = "rawdata/{sample}"
    log: "logs/{sample}/fastqdump.log"
    output:
        r1 = temp("rawdata/{sample}_1.fastq"),
        r2 = temp("rawdata/{sample}_2.fastq")
    conda:
        "envs/deeptools.yml"
    shell: "fastq-dump --split-files --outdir {params.outdir} {input} 2> {log}"

rule TrimReads:
    input:
        r1 = "rawdata/{sample}_R1_001.fastq.gz",
        r2 = "rawdata/{sample}_R2_001.fastq.gz"
    params:
        trimmomaticdir = config["trimmomaticdir"]
    threads: 2
    log: "logs/{sample}/TrimReads.log"
    output:
        r1 = "rawdata/{sample}_1.trimmed.fastq.gz",
        r2 = "rawdata/{sample}_2.trimmed.fastq.gz",
        utr1 = "rawdata/{sample}_1.trimmed.fastq.unpaired.gz",
        utr2 = "rawdata/{sample}_2.trimmed.fastq.unpaired.gz"
    conda:
        "envs/deeptools.yml"
    shell:
        """
java -jar {params.trimmomaticdir}/trimmomatic-0.39.jar \
PE -phred33 {input.r1} \
{input.r2} \
{output.r1} \
{output.utr1} \
{output.r2} \
{output.utr2} \
ILLUMINACLIP:{params.trimmomaticdir}/adapters/TruSeq3-PE.fa:2:30:10:8:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
2> {log}
        """
rule AlignBowtie2:
    input:
        r1 = "rawdata/{sample}_1.trimmed.fastq.gz",
        r2 = "rawdata/{sample}_2.trimmed.fastq.gz"
    params:
        index = config["bowtie2hs1index"]
    output:
        sam = temp("outs/{sample}/{sample}.sam")
    log:
        out = "logs/{sample}/bowtie2.out",
        err = "logs/{sample}/bowtie2.err"
    threads: 2
    conda:
        "envs/deeptools.yml"
    shell:
        "bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.sam} --threads {threads} > {log.out} 2> {log.err}"

rule alignSTAR:
    input:
        r1 = "rawdata/{sample}_1.trimmed.fastq.gz",
        r2 = "rawdata/{sample}_2.trimmed.fastq.gz"
    params:
        index = config["starindex"],
        outdir = "outs/{sample}/star_output/"
    output:
        bam = temp("outs/{sample}/star_output/Aligned.out.sam")
    log:
        out = "logs/{sample}/STAR.out",
        err = "logs/{sample}/STAR.err"
    threads: 10
    shell:
        """
STAR --genomeDir {params.index} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {params.outdir} --runThreadN {threads} --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 > {log.out} 2> {log.err}
        """

rule sortSTARbams:
    input:
        bamSTAR = "outs/{sample}/star_output/Aligned.out.sam"
    output:
        sortedbamSTAR = "outs/{sample}/star_output/{sample}.STAR.sorted.bam"
    log: "logs/{sample}/sortSTARbams.log"
    threads: 6
    conda:
        "envs/deeptools.yml"
    shell: "samtools sort {input.bamSTAR} -o {output.sortedbamSTAR} 2> {log}"

rule indexsortedSTARbam:
    input:
        sortedbamSTAR = "outs/{sample}/star_output/{sample}.STAR.sorted.bam"
    log: "logs/{sample}/indexsortedSTARbam.log"
    output:
        sortedbamSTARindex = "outs/{sample}/star_output/{sample}.STAR.sorted.bam.bai"
    threads: 6
    conda:
        "envs/deeptools.yml"
    shell: "samtools index {input.sortedbamSTAR} 2> {log}"


rule samtobam:
    input:
        sam = "outs/{sample}/{sample}.sam"
    output:
        bam = pipe("outs/{sample}/{sample}.bam")
    threads: 2
    conda:
        "envs/deeptools.yml"
    shell:
        "samtools view -S -b {input.sam} > {output.bam}"

rule sortbams:
    input:
        bam = "outs/{sample}/{sample}.bam"
    output:
        sortedbam = "outs/{sample}/{sample}.sorted.bam"
    threads: 2
    conda:
        "envs/deeptools.yml"
    shell: "samtools sort {input.bam} -o {output.sortedbam}"

rule indexsortedbam:
    input:
        sortedbam = "outs/{sample}/{sample}.sorted.bam"
    log: "logs/{sample}/indexsortedbam.log"
    output:
        bamindex = "outs/{sample}/{sample}.sorted.bam.bai"
    threads: 2
    conda:
        "envs/deeptools.yml"
    shell: "samtools index {input.sortedbam} 2> {log}"

rule bamstats:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam"
    output:
        bamstats = "outs/{sample}/{sample}.bam.sorted.stats.txt"
    threads: 2
    conda:
        "envs/deeptools.yml"
    shell:
        "samtools stats {input.bam} > {output.bamstats}"

#######################################################################
#t2t style analysis of repeats
#######################################################################
rule repeatCounts:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam",
        bamindex = "outs/{sample}/{sample}.sorted.bam.bai",
        repeats = config["repeats"]
    log: "logs/{sample}/repeatCounts.log"
    output:
        repeatcounts = "outs/{sample}/{sample}.repeatmasker.counts"
    threads: 2
    conda:
        "envs/deeptools.yml"
    shell:
        "bedtools coverage -counts -F 0.5 -sorted -b {input.bam} -a {input.repeats} > {output.repeatcounts} 2> {log}"


rule countMappednonmitoreads:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam"
    output:
        mnma = "outs/{sample}/{sample}_number_nonmito_mapped_per_million.txt"
    log: "logs/{sample}/countMappednonmitoreads.log"
    shell:
        "samtools view -F 4 /home/mk/scienceL/marco/outs/SRR6515349/SRR6515349.sorted.bam | awk '$3 !~ /^MT/ { count++ } END { print count/1000000 }' > {output.mnma} "

rule getnormalizedCounts:
    input:
        repeatcounts = "outs/{sample}/{sample}.repeatmasker.counts",
        mnma = "outs/{sample}/{sample}_number_nonmito_mapped_per_million.txt"
    threads: 2
    output:
        repeatcounts = "results/{sample}/{sample}.repeatmasker.countsPerMillionNonMito"
    log: "logs/{sample}/getnormalizedCounts.log"
    conda:
        "envs/deeptools.yml"
    shell:
        """
        awk -v milnonmito=$(cat {input.mnma}) -F "\t" '{{print $0 "\t" $NF/milnonmito}}' {input.repeatcounts} > {output.repeatcounts} 2> {log}
        """

rule getnormalizedFamilyCounts:
    input:
        counts = "results/{sample}/{sample}.repeatmasker.countsPerMillionNonMito"
    threads: 2
    output:
        l1 = "results/{sample}/{sample}.l1.counts.txt",
        l1hs = "results/{sample}/{sample}.l1hs.counts.txt",
        alu = "results/{sample}/{sample}.alu.counts.txt",
        aluy = "results/{sample}/{sample}.aluy.counts.txt",
        herv = "results/{sample}/{sample}.herv.counts.txt",
        hervk = "results/{sample}/{sample}.hervk.counts.txt",
        sva = "results/{sample}/{sample}.sva.counts.txt"
    conda:
        "envs/deeptools.yml"
    shell:
        """
        awk -F "\t" '/LINE\/L1/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.l1}
awk -F "\t" '/L1HS/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.l1hs}
awk -F "\t" '/SINE\/Alu/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.alu}
awk -F "\t" '/AluY/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.aluy}
awk -F "\t" '/LTR\/ERV/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.herv}
awk -F "\t" '/HERVK/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.hervk}
awk -F "\t" '/SVA/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.sva}
        """

rule multiBamSummary:
    input:
        bam = expand("outs/{a}/{a}.sorted.bam", a=samples)
    output:
        agg = "results/agg/multiBamSummary.npz"
    threads: 2
    log: "logs/multiBamSummary.log"
    conda:
        "envs/deeptools.yml"
    shell: "multiBamSummary bins --numberOfProcessors {threads} --bamfiles {input.bam} -o {output.agg} 2> {log}"

rule plotbamPCA:
    input:
        agg = "results/agg/multiBamSummary.npz"
    output:
        plot = report("results/agg/plots/PCA_readCounts.png", caption = "report/plotbamPCA.rst", category="PCA")
    log: "logs/agg/plotbamPCA.log"
    conda:
        "envs/deeptools.yml"
    shell:
        "plotPCA -in {input.agg} -o {output.plot} -T 'PCA of read counts' 2> {log}"

#######################################################################
#Get counts via several methods
#######################################################################
rule TElocal:
    input:
        sortedSTARbam = "outs/{sample}/star_output/{sample}.STAR.sorted.bam",
        sortedbamSTARindex = "outs/{sample}/star_output/{sample}.STAR.sorted.bam.bai"
    log: "logs/{sample}/TElocal.log"
    params:
        refseq = config["refseq"],
        locindTElocal = config["locindTElocal"],
        outputprefixMULTI = "outs/{sample}/TElocal/{sample}_MULTI",
        outputprefixUNIQUE = "outs/{sample}/TElocal/{sample}_UNIQUE"
    output:
        countsMULTI = "outs/{sample}/TElocal/{sample}_MULTI.cntTable",
        countsUNIQUE = "outs/{sample}/TElocal/{sample}_UNIQUE.cntTable"
    threads: 7
    shell: 
        """
TElocal --sortByPos -b {input.sortedSTARbam} --GTF {params.refseq} --TE {params.locindTElocal} --project {params.outputprefixMULTI} 2> {log}
TElocal --sortByPos -b {input.sortedSTARbam} --GTF {params.refseq} --TE {params.locindTElocal} --mode uniq --project {params.outputprefixUNIQUE} 2>> {log}
        """


rule featureCounts:
    input:
        sortedSTARbams = expand("outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = samples)
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
#######################################################################
# deeptools plots
#######################################################################

rule deeptools_coverage:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam"
    conda:
        "envs/deeptools.yml"
    output:
        coverage = "outs/{sample}/{sample}_cov.bw"
    shell: "bamCoverage -b {input.bam} -o {output.coverage} --outFileFormat bigwig"

rule deeptools_plotcoverage:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam"
    params:
        l1hs6kb = config["l1hs6kb"],
        l1hs6kbintact = config["l1hs6kbintact"],
    conda:
        "envs/deeptools.yml"
    output:
        genomecoverageplot = "results/{sample}/plots/coverageGenome",
        l1hscoverageplot = "results/{sample}/plots/coverageL1HS",
        l1hsintactcoverageplot = "results/{sample}/plots/coverageL1HSintact"
    shell:
        """
plotCoverage -b {input.bam} \
--plotFile {output.genomecoverageplot} \
-n 1000000 \
--plotTitle "Whole Genome Coverage" \
--ignoreDuplicates \
--minMappingQuality 10

plotCoverage -b {input.bam} \
--BED {params.l1hs6kb} \
--plotFile {output.l1hscoverageplot} \
-n 1000000 \
--plotTitle "L1HS Coverage" \
--ignoreDuplicates \
--minMappingQuality 10

plotCoverage -b {input.bam} \
--BED {params.l1hs6kbintact} \
--plotFile {output.l1hsintactcoverageplot} \
-n 1000000 \
--plotTitle "L1HS Coverage" \
--ignoreDuplicates \
--minMappingQuality 10
        """



rule deeptools_plotAggcoverage:
    input:
        coverage = "outs/{sample}/coverage.bw"
    params:
        l1hs6kb = config["l1hs6kb"],
        l1hs6kbintact = config["l1hs6kbintact"]
    conda:
        "envs/deeptools.yml"
    log: "logs/{sample}/deeptools_plotAggcoverage.log"
    output:
        matrix = "outs/{sample}/matrix_l1hs6kb_coverage.tab.gz",
        heatmap_coverage_k3 = report("results/{sample}/plots/heatmap_coverage_k3.png", caption = "report/deeptools_plotAggcoverageheatmap_coverage_k3.rst", category="coverage"),
        heatmap_coverage = report("results/{sample}/plots/heatmap_coverage.png", caption = "report/deeptools_plotAggcoverageheatmap_coverage.rst", category="coverage"),
        matrixintact = "outs/{sample}/matrix_l1hs6kbintact_coverage.tab.gz",
        heatmap_coverage_k3intact = report("results/{sample}/plots/heatmap_coverageintact_k3.png", caption = "report/deeptools_plotAggcoverageheatmap_coverageintact_k3.rst", category="coverage"),
        heatmap_coverageintact = report("results/{sample}/plots/heatmap_coverageintact.png", caption = "report/deeptools_plotAggcoverageheatmap_coverageintact.rst", category="coverage")

    shell:
        """
computeMatrix reference-point \
 -S {input.coverage} \
 -R {params.l1hs6kb} \
 --referencePoint TSS \
 --upstream 3000 \
 --binSize 100 \
 --downstream 9000 \
 -out {output.matrix} 2> {log}

 plotHeatmap \
 -m {output.matrix} \
 -out {output.heatmap_coverage_k3} \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --sortRegions descend \
 --sortUsing sum \
 --missingDataColor 0 \
 --labelRotation 0 \
 --linesAtTickMarks \
 --kmeans 3 \
 --samplesLabel "L1Hs6kb Coverage" \
 --dpi 400 \
 --colorMap Blues \
 --plotTitle "" 2>> {log}

 plotHeatmap \
 -m {output.matrix} \
 -out {output.heatmap_coverage} \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --sortRegions descend \
 --sortUsing sum \
 --regionsLabel l1hs \
 --missingDataColor 0 \
 --labelRotation 0 \
 --linesAtTickMarks \
 --samplesLabel "L1Hs6kb Coverage" \
 --dpi 400 \
 --colorMap Blues \
 --plotTitle "" 2>> {log}

 computeMatrix reference-point \
 -S {input.coverage} \
 -R {params.l1hs6kbintact} \
 --referencePoint TSS \
 --upstream 3000 \
 --binSize 100 \
 --downstream 9000 \
 -out {output.matrixintact} 2>> {log}

 plotHeatmap \
 -m {output.matrixintact} \
 -out {output.heatmap_coverage_k3intact} \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --sortRegions descend \
 --sortUsing sum \
 --missingDataColor 0 \
 --labelRotation 0 \
 --linesAtTickMarks \
 --kmeans 3 \
 --samplesLabel "L1Hs6kb intact Coverage" \
 --dpi 400 \
 --colorMap Blues \
 --plotTitle "" 2>> {log}

 plotHeatmap \
 -m {output.matrixintact} \
 -out {output.heatmap_coverageintact} \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --sortRegions descend \
 --sortUsing sum \
 --regionsLabel l1hs \
 --missingDataColor 0 \
 --labelRotation 0 \
 --linesAtTickMarks \
 --samplesLabel "L1Hs6kb intact Coverage" \
 --dpi 400 \
 --colorMap Blues \
 --plotTitle "" 2>> {log}
        """

#######################################################################
#DE analysis
#######################################################################
rule mergeTElocal:
    input:
        telocalMULTI = expand("outs/{sample}/TElocal/{sample}_MULTI.cntTable", sample = samples),
        telocalUNIQUE = expand("outs/{sample}/TElocal/{sample}_UNIQUE.cntTable", sample = samples)

    params:
        sampleinfo = config["sample_table"]
    conda: "envs/renv.yml"
    log: "logs/agg/mergeTElocal.log"
    output:
        aggcountsMULTI = "outs/agg/TElocalCounts_MULTI.txt",
        aggcountsUNIQUE = "outs/agg/TElocalCounts_UNIQUE.txt"
    script:
        "scripts/mergeTEdf.R"


rule DEseq2:
    input:
        star = "outs/agg/refseq.counts.txt",
        telocal_MULTI = "outs/agg/TElocalCounts_MULTI.txt",
        telocal_UNIQUE = "outs/agg/TElocalCounts_UNIQUE.txt",
    params:
        sample_table = config["sample_table"],
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        levels = config["levels"],
        outputdir = "results/agg/deseq2"
    conda: "envs/renv.yml"
    log: "logs/agg/DEseq2.log"
    output:
        results = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"])
    script:
        "scripts/DESeq2.R"

rule clusterprofiler:
    input:
        deresults = expand("results/agg/deseq2/star/{contrast}/results.csv", contrast = config["contrasts"])
    params:
        contrasts = config["contrasts"],
        SenMayoHuman = config["SenMayoHuman"],
        inputdir = "results/agg/deseq2/star",
        outputdir = "results/agg/clusterprofiler"
    conda:
        "envs/Rclusterprofiler.yml"
    log:
        "logs/agg/clusterprofiler.log"
    output:
        enrichplotsUP = report(expand("results/agg/clusterprofiler/{contrast}/hypgeo/go_enrichedUP.pdf", contrast = config["contrasts"]), caption = "report/clusterprofilergo_enrichedUP.rst", category="enrichment"),
        enrichplotsDOWN = report(expand("results/agg/clusterprofiler/{contrast}/hypgeo/go_enrichedDOWN.pdf", contrast = config["contrasts"]), caption = "report/clusterprofilergo_enrichedDOWN.rst", category="enrichment"),
        SenMayoPlot = report(expand("results/agg/clusterprofiler/{contrast}/gsea/SenMayo.pdf", contrast = config["contrasts"]),caption = "report/clusterprofilerSenMayo.rst", category="enrichment"),
        ridge = report(expand("results/agg/clusterprofiler/{contrast}/gsea/ridgeplot.pdf", contrast = config["contrasts"]),caption = "report/clusterprofilerridgeplot.rst", category="enrichment"),
        dotplot = report(expand("results/agg/clusterprofiler/{contrast}/gsea/dotplot.pdf", contrast = config["contrasts"]),caption = "report/clusterprofilerdotplot.rst", category="enrichment"),
        emapplot = report(expand("results/agg/clusterprofiler/{contrast}/gsea/emapplot.pdf", contrast = config["contrasts"]),caption = "report/clusterprofileremapplot.rst", category="enrichment"),
        outfile = "results/agg/clusterprofiler/outfile.txt"
    script:
        "scripts/clusterprofiler.R"



rule pyGenomeTracks:
    input:
        bw = expand("outs/{sample}/{sample}_cov.bw", sample = samples)
    params:
        refseq = config["refseq"],
        l1hs6kbintactbed = config["l1hs6kbintactbed"],
        repeatsbed = config["repeatsbed"],
        outputdir = "results/agg/repeatanalysis/genometracks"
    conda:
        "envs/deeptools.yml"
    log:
        "logs/agg/pyGenomeTracks.log"
    output:
        outfile = "results/agg/repeatanalysis/genometracks/outfile.txt"
    shell:
        """
make_tracks_file --trackFiles \
{params.l1hs6kbintactbed} \
{input.bw} \
{params.repeatsbed} \
-o {params.outputdir}/atracks.ini 2> {log}

#modify tracks ini file to show labels
sed 's/labels = false/labels = true/g' {params.outputdir}/atracks.ini > {params.outputdir}/atracksMOD.ini

#manually modify tracks to show labels
cat {params.l1hs6kbintactbed} | while read line
do
chr=$(awk '{{print $1}}' <<< $line)
start=$(awk '{{print $2}}' <<< $line)
stop=$(awk '{{print $3}}' <<< $line)
echo $chr $start $stop
pyGenomeTracks --tracks {params.outputdir}/atracksMOD.ini --region ${{chr}}:$((${{start}}-500))-$((${{stop}}+500)) \
--dpi 300 -o {params.outputdir}/${{chr}}${{start}}${{stop}}.png 2>> {log}
done

touch {output.outfile}
        """



rule repeatanalysis:
    input:
        expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"])
    params:
        contrasts = config["contrasts"],
        telocaltypes = config["telocaltypes"],
        sample_table = config["sample_table"],
        repeats = config["repeats"],
        inputdir = "results/agg/deseq2",
        outputdir = "results/agg/repeatanalysis"
    conda:
        "envs/repeatanalysis.yml"
    log:
        "logs/agg/repeatanalysis.log"
    output:
        activeelementContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/activeelementContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysisactiveelementContrastplot.rst", category="repeat analysis"),
        familyContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/familyContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysisfamilyContrastplot.rst", category="repeat analysis"),
        combinedelementContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/combinedContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysiscombinedContrastplot.rst", category="repeat analysis"),
        outfile = "results/agg/repeatanalysis/outfile.txt"
    script:
        "scripts/repeatanalysis.R"

rule deeptools_plotAggSignal:
    input:
        coverage = expand("outs/{sample}/{sample}_cov.bw", sample = samples)
    params:
        l1hs6kb = config["l1hs6kb"],
        aluy = config["aluY"],
        hervk = config["HERVK"],
        l1hs6kbintact = config["l1hs6kbintact"],
        outputdirprofile = "results/agg/repeatanalysis",
        outputdirmatrix = "outs/agg"
    conda:
        "envs/deeptools.yml"
    log: "logs/agg/deeptools_plotAggSignal.log"
    output:
        matrixl1hs = "outs/agg/matrix_l1hs6kb_coverage.tab.gz",
        profilel1hs = "results/agg/repeatanalysis/L1hs_profile.png",
        profileheatl1hs = "results/agg/repeatanalysis/L1hs_profile_heatmap.png",
        matrixaluy = "outs/agg/matrix_aluy_coverage.tab.gz",
        profilealuy = "results/agg/repeatanalysis/aluy_profile.png",
        profileheataluy = "results/agg/repeatanalysis/aluy_profile_heatmap.png",
        matrixhervk = "outs/agg/matrix_hervk_coverage.tab.gz",
        profilehervk = "results/agg/repeatanalysis/hervk_profile.png",
        profileheathervk = "results/agg/repeatanalysis/hervk_profile_heatmap.png",
    shell:
        """
computeMatrix reference-point \
 -S {input.coverage} \
 -R {params.l1hs6kb} \
 --referencePoint TSS \
 --upstream 3000 \
 --binSize 100 \
 --downstream 9000 \
 -out {output.matrixl1hs} 2> {log}

 plotProfile \
 -m {output.matrixl1hs} \
 --perGroup \
 -out {output.profilel1hs} \
 --dpi 400 \
 --plotTitle "L1HS" 2>> {log}

  plotProfile \
 -m {output.matrixl1hs} \
 --perGroup \
 --plotType heatmap \
 -out {output.profileheatl1hs} \
 --dpi 400 \
 --plotTitle "L1HS" 2>> {log}

 ##
 computeMatrix reference-point \
 -S {input.coverage} \
 -R {params.aluy} \
 --referencePoint TSS \
 --upstream 3000 \
 --binSize 100 \
 --downstream 9000 \
 -out {output.matrixaluy} 2> {log}

 plotProfile \
 -m {output.matrixaluy} \
 --perGroup \
 -out {output.profilealuy} \
 --dpi 400 \
 --plotTitle "AluY" 2>> {log}

  plotProfile \
 -m {output.matrixaluy} \
 --perGroup \
 --plotType heatmap \
 -out {output.profileheataluy} \
 --dpi 400 \
 --plotTitle "AluY" 2>> {log}
 ##
 computeMatrix reference-point \
 -S {input.coverage} \
 -R {params.hervk} \
 --referencePoint TSS \
 --upstream 3000 \
 --binSize 100 \
 --downstream 9000 \
 -out {output.matrixhervk} 2> {log}

 plotProfile \
 -m {output.matrixhervk} \
 --perGroup \
 -out {output.profilehervk} \
 --dpi 400 \
 --plotTitle "HERVK" 2>> {log}

  plotProfile \
 -m {output.matrixhervk} \
 --perGroup \
 --plotType heatmap \
 -out {output.profileheathervk} \
 --dpi 400 \
 --plotTitle "HERVK" 2>> {log}
    """



#  --heatmapHeight 15  \
#  --refPointLabel TSS \
#  --sortRegions descend \
#  --sortUsing sum \
#  --missingDataColor 0 \
#  --labelRotation 0 \
#  --linesAtTickMarks \
#  --kmeans 3 \
#  --samplesLabel "L1Hs6kb Coverage" \

