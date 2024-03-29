import os
import pandas as pd
from pathlib import Path
configfile: "conf/private/configPrivate.yaml"
configfile: "conf/shared/configShared.yaml"
pepfile: "conf/private/project_config.yaml"



container: "docker://maxfieldkelsey/snakesenescence:latest"

samples = pep.sample_table.sample_name
samples = ["PRO1", "PRO2", "PRO3", "NT1", "NT2", "NT3"]
peptable = pep.sample_table
import csv
peptable.to_csv("conf/private/peptable.csv", index = False, quoting=csv.QUOTE_NONNUMERIC)

telocaltypes = config["telocaltypes"]
contrasts = config["contrasts"]
counttypes = config["counttypes"]
directions = ["UP", "DOWN"]
def exHelper(path):
    try:
        files = expand(path, telocaltype = telocaltypes, contrast = contrasts, rte = ["L1HS", "AluY", "HERVK-int"], direction = ["UP", "DOWN"], sample = samples)
        return files
    except:
        None
    try:
        files = expand(path, telocaltype = telocaltypes, contrast = contrasts, direction = ["UP", "DOWN"])
        return files
    except:
        None
    try:
        files = expand(path, telocaltype = telocaltypes, contrast = contrasts, direction = ["UP", "DOWN"])
        return files
    except:
        None
    try:
        files = expand(path, telocaltype = telocaltypes, contrast = contrasts)
        return files
    except:
        None
    try:
        files = expand(path, telocaltype = telocaltypes)
        return files
    except:
        None

# tips:
# build needed directories in python! Will save you much frustration :)
# watch out for multiline shell commands - they fail for weird reasons sometimes.
try:
    dertes = pd.read_table("results/agg/repeatanalysis/allactiveDETEs.tsv", header = None)
except:
    None
try:
    gvizout = expand("temp/outfile{telocaltype}{contrast}{rtekind}{range}._{telocaltype}._{contrast}._{rtekind}._{range}.gviz.txt", 
        telocaltype = config["telocaltypes"],
        contrast = config["contrasts"],
        rtekind = list(dertes.iloc[:,0].unique()),
        range = [500, 1000, 6000, 30000])
except:
    gvizout = ""

rule all:
    input:
        expand("outs/{sample}/star_output/{sample}_forward_BPM_cov.bw", sample = samples)
        # expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.alignment.pdf", telocaltype = telocaltypes, contrast = contrasts, rte = ["L1HS", "AluY", "HERVK-int"], direction = ["UP", "DOWN"])
        
########################################################### Make folders
if True:
    path = 'slurm'
    if not os.path.exists(path):
        os.makedirs(path)

    path = 'templates'
    if not os.path.exists(path):
        os.makedirs(path)

    path = 'tools'
    if not os.path.exists(path):
        os.makedirs(path)

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

    path = "results/agg/repeatanalysis/genometracks"
    if not os.path.exists(path):
        os.makedirs(path)


    for counttype in config["counttypes"]:
        path = "results/agg/repeatanalysis/%s"%(counttype)
        if not os.path.exists(path):
            os.makedirs(path)

    for counttype in config["counttypes"]:
        for contrast in config["contrasts"]:
            path = 'results/agg/repeatanalysis/%s/%s'%(counttype,contrast)
            if not os.path.exists(path):
                os.makedirs(path)

    for counttype in config["counttypes"]:
        path = "results/agg/repeatanalysis/genometracks/%s"%(counttype)
        if not os.path.exists(path):
            os.makedirs(path)

    for counttype in config["counttypes"]:
        for contrast in config["contrasts"]:
            path = 'results/agg/repeatanalysis/genometracks/%s/%s'%(counttype,contrast)
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

    if not os.path.exists('results/agg/clusterprofiler/genes'):
        os.makedirs('results/agg/clusterprofiler/genes')

    for contrast in config["contrasts"]:
        if not os.path.exists('results/agg/clusterprofiler/%s'%(contrast)):
            os.makedirs('results/agg/clusterprofiler/%s'%(contrast))

    for contrast in config["contrasts"]:
        if not os.path.exists('results/agg/clusterprofiler/%s/genes'%(contrast)):
            os.makedirs('results/agg/clusterprofiler/%s/genes'%(contrast))

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

        for contrast in config["contrasts"]:
            path = 'results/agg/clusterprofiler/%s/KEGG/topsigpathways'%(contrast)
            if not os.path.exists(path):
                os.makedirs(path)

        for contrast in config["contrasts"]:
            path = 'results/agg/clusterprofiler/%s/KEGG/importantpathways'%(contrast)
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

        outpath = 'outs/%s/trimmedReads' % (e)
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
    threads: 4
    conda:
        "deeptools"
    shell: "prefetch {wildcards.sample} --output-directory rawdata"

rule fastqdump:
    input: "rawdata/{sample}/{sample}.sra"
    threads: 4
    params:
        outdir = "rawdata/{sample}"
    log: "logs/{sample}/fastqdump.log"
    output:
        r1 = temp("rawdata/{sample}_1.fastq"),
        r2 = temp("rawdata/{sample}_2.fastq")
    conda:
        "deeptools"
    shell: "fastq-dump --split-files --outdir {params.outdir} {input} 2> {log}"

        
# rule TrimReads:
#     input:
#         r1=lambda wildcards: peptable[peptable["sample_name"] == '{sample}'.format(sample=wildcards.sample)]["file_path_R1"],
#         r2=lambda wildcards: peptable[peptable["sample_name"] == '{sample}'.format(sample=wildcards.sample)]["file_path_R2"]
#         # r1 = "rawdata/{sample}_1.fastq.gz",
#         # r2 = "rawdata/{sample}_2.fastq.gz"
#     params:
#         trimmomaticdir = config["trimmomaticdir"]
#     threads: 2
#     log: "logs/{sample}/TrimReads.log"
#     output:
#         r1 = "rawdata/{sample}_1.trimmed.fastq.gz",
#         r2 = "rawdata/{sample}_2.trimmed.fastq.gz",
#         utr1 = "rawdata/{sample}_1.trimmed.fastq.unpaired.gz",
#         utr2 = "rawdata/{sample}_2.trimmed.fastq.unpaired.gz"
#     conda:
#         "deeptools"
#     shell:
#         """
# java -jar {params.trimmomaticdir}/trimmomatic-0.39.jar \
# PE -phred33 {input.r1} \
# {input.r2} \
# {output.r1} \
# {output.utr1} \
# {output.r2} \
# {output.utr2} \
# ILLUMINACLIP:{params.trimmomaticdir}/adapters/TruSeq3-PE.fa:2:30:10:8:true \
# LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
# 2> {log}
#         """


rule fastp:
    input:
        r1=lambda wildcards: peptable[peptable["sample_name"] == '{sample}'.format(sample=wildcards.sample)]["file_path_R1"],
        r2=lambda wildcards: peptable[peptable["sample_name"] == '{sample}'.format(sample=wildcards.sample)]["file_path_R2"]
    params:
    threads: 6
    log: "logs/{sample}/fastp.log"
    output:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
        r2 = "outs/{sample}/trimmedReads/{sample}_2.trimmed.fastq.gz",
        json = "outs/{sample}/trimmedReads/fastp.json",
        html = "outs/{sample}/trimmedReads/fastp.html"
    conda:
        "qc"
    shell:
        """
fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --json {output.json} --html {output.html} --detect_adapter_for_pe --thread {threads}
        """

#"rawdata/X3TC_2.trimmed.fastq.gz"
###################### Bowtie2
rule AlignBowtie2:
    input:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
        r2 = "outs/{sample}/trimmedReads/{sample}_2.trimmed.fastq.gz"
    params:
        index = config["bowtie2hs1index"]
    output:
        sam = temp("outs/{sample}/{sample}.sam")
    log:
        out = "logs/{sample}/bowtie2.out",
        err = "logs/{sample}/bowtie2.err"
    threads: 4
    conda:
        "deeptools"
    shell:
        "bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.sam} --threads {threads} > {log.out} 2> {log.err}"

rule samtobam:
    input:
        sam = "outs/{sample}/{sample}.sam"
    output:
        bam = temp("outs/{sample}/{sample}.bam")
    threads: 4
    conda:
        "deeptools"
    shell:
        "samtools view -S -b {input.sam} > {output.bam}"

rule sortbams:
    input:
        bam = "outs/{sample}/{sample}.bam"
    output:
        sortedbam = "outs/{sample}/{sample}.sorted.bam"
    threads: 4
    conda:
        "deeptools"
    shell: "samtools sort {input.bam} -o {output.sortedbam}"

rule indexsortedbam:
    input:
        sortedbam = "outs/{sample}/{sample}.sorted.bam"
    log: "logs/{sample}/indexsortedbam.log"
    output:
        bamindex = "outs/{sample}/{sample}.sorted.bam.bai"
    threads: 4
    conda:
        "deeptools"
    shell: "samtools index {input.sortedbam} 2> {log}"

rule bamstats:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam"
    output:
        bamstats = "outs/{sample}/{sample}.bam.sorted.stats.txt"
    threads: 4
    conda:
        "deeptools"
    shell:
        "samtools stats {input.bam} > {output.bamstats}"

############ STAR


rule alignSTAR:
    input:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
        r2 = "outs/{sample}/trimmedReads/{sample}_2.trimmed.fastq.gz"
    params:
        index = config["starindex"],
        outdir = "outs/{sample}/star_output/"
    output:
        bam = temp("outs/{sample}/star_output/Aligned.out.sam")
    log:
        out = "logs/{sample}/STAR.out",
        err = "logs/{sample}/STAR.err"
    threads: 8
    resources:
        mem_mb  = 60000
    conda:
        "deeptools"
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
    threads: 4
    conda:
        "deeptools"
    shell: "samtools sort {input.bamSTAR} -o {output.sortedbamSTAR} 2> {log}"

rule indexsortedSTARbam:
    input:
        sortedbamSTAR = "outs/{sample}/star_output/{sample}.STAR.sorted.bam"
    log: "logs/{sample}/indexsortedSTARbam.log"
    output:
        sortedbamSTARindex = "outs/{sample}/star_output/{sample}.STAR.sorted.bam.bai"
    threads: 4
    conda:
        "deeptools"
    shell: "samtools index {input.sortedbamSTAR} 2> {log}"


#######################################################################
#Bowtie based alignment stats
#######################################################################
rule repeatCounts:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam",
        bamindex = "outs/{sample}/{sample}.sorted.bam.bai",
        repeats = config["repeats"]
    log: "logs/{sample}/repeatCounts.log"
    output:
        repeatcounts = "outs/{sample}/{sample}.repeatmasker.counts"
    threads: 4
    conda:
        "deeptools"
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
    threads: 4
    output:
        repeatcounts = "results/{sample}/{sample}.repeatmasker.countsPerMillionNonMito"
    log: "logs/{sample}/getnormalizedCounts.log"
    conda:
        "deeptools"
    shell:
        """
        awk -v milnonmito=$(cat {input.mnma}) -F "\t" '{{print $0 "\t" $NF/milnonmito}}' {input.repeatcounts} > {output.repeatcounts} 2> {log}
        """

rule getnormalizedFamilyCounts:
    input:
        counts = "results/{sample}/{sample}.repeatmasker.countsPerMillionNonMito"
    threads: 4
    output:
        l1 = "results/{sample}/{sample}.l1.counts.txt",
        l1hs = "results/{sample}/{sample}.l1hs.counts.txt",
        alu = "results/{sample}/{sample}.alu.counts.txt",
        aluy = "results/{sample}/{sample}.aluy.counts.txt",
        herv = "results/{sample}/{sample}.herv.counts.txt",
        hervk = "results/{sample}/{sample}.hervk.counts.txt",
        sva = "results/{sample}/{sample}.sva.counts.txt"
    conda:
        "deeptools"
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
    threads: 4
    log: "logs/multiBamSummary.log"
    conda:
        "deeptools"
    shell: "multiBamSummary bins --numberOfProcessors {threads} --bamfiles {input.bam} -o {output.agg} 2> {log}"

rule plotbamPCA:
    input:
        agg = "results/agg/multiBamSummary.npz"
    output:
        plot = report("results/agg/plots/PCA_readCounts.png", caption = "report/plotbamPCA.rst", category="PCA")
    log: "logs/agg/plotbamPCA.log"
    conda:
        "deeptools"
    shell:
        "plotPCA -in {input.agg} -o {output.plot} -T 'PCA of read counts' 2> {log}"

#######################################################################
#Get counts via several methods
#######################################################################
#note that maptype must be either uniq or multi
rule TElocal:
    input:
        sortedSTARbam = "outs/{sample}/star_output/{sample}.STAR.sorted.bam",
        sortedbamSTARindex = "outs/{sample}/star_output/{sample}.STAR.sorted.bam.bai"
    log: "logs/{sample}/TElocal_{maptype}.log"
    params:
        refseq = config["refseq"],
        locindTElocal = config["locindTElocal"],
        telocalstrandparam = config['telocalstrandparam'],
        outputprefix = "outs/{sample}/TElocal/{sample}_{maptype}",
    output:
        counts = "outs/{sample}/TElocal/{sample}_{maptype}.cntTable",
    threads: 4
    resources:
        mem_mb  = 40000
    shell: 
        """
TElocal --sortByPos -b {input.sortedSTARbam} --stranded {params.telocalstrandparam} --mode {wildcards.maptype} --GTF {params.refseq} --TE {params.locindTElocal} --project {params.outputprefix} 2> {log}
        """


rule featureCounts:
    input:
        sortedSTARbams = expand("outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = samples)
    output:
        countsmessy = "outs/agg/refseq.counts_messy.txt",
        counts = "outs/agg/refseq.counts.txt",
        countsstrandnonspecificmessy = "outs/agg/refseq.countsstrandnonspecific_messy.txt",
        countsstrandnonspecific = "outs/agg/refseq.countsstrandnonspecific.txt",
        metafeaturecounts = "outs/agg/refseq.metafeature.counts.txt"
    params: 
        gtf = config['refseq'],
        featureCountsstrandparam = config['featureCountsstrandparam']
    log: "logs/agg/featureCounts.log"
    conda:
        "deeptools"
    threads: 4
    shell: 
        """
featureCounts -p -s {params.featureCountsstrandparam} -T {threads} -t exon -a {params.gtf} -o {output.countsmessy} {input.sortedSTARbams} 2> {log}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
featureCounts -p -s {params.featureCountsstrandparam} -T {threads} -B -O -a {params.gtf} -o {output.countsstrandnonspecificmessy} {input.sortedSTARbams} 2>> {log}
cut -f1,7- {output.countsstrandnonspecificmessy} | awk 'NR > 1' > {output.countsstrandnonspecific}
featureCounts -p -T {threads} -B -O -a {params.gtf} -o {output.metafeaturecounts} {input.sortedSTARbams} 2>> {log}
        """
#######################################################################
# deeptools plots
#######################################################################

rule normalized_stranded_coverage:
    input:
        bam = "outs/{sample}/star_output/{sample}.STAR.sorted.bam",
        bamindex = "outs/{sample}/star_output/{sample}.STAR.sorted.bam.bai"
    params:
        hs1size = config["hs1size"],
        effectivegenomesize = 3117292070 #derived from using faCount from kenttools and summing all non N bases in genome
    conda:
        "deeptools"
    resources:
        cpus_per_task = 10
    output:
        bwF = "outs/{sample}/star_output/{sample}_forward_BPM_cov.bw",
        bwR = "outs/{sample}/star_output/{sample}_reverse_BPM_cov.bw"
    log: "logs/{sample}/normalized_stranded_coverage.log"
    shell:
        """
#forward
bamCoverage \
-b {input.bam} \
-o {output.bwF} \
-of bigwig \
--filterRNAstrand forward \
--numberOfProcessors max \
--normalizeUsing BPM \
--binSize 10 \
--effectiveGenomeSize {params.effectivegenomesize}

#reverse
bamCoverage \
-b {input.bam} \
-o {output.bwR} \
-of bigwig \
--filterRNAstrand reverse \
--numberOfProcessors max \
--normalizeUsing BPM \
--binSize 10 \
--effectiveGenomeSize {params.effectivegenomesize}
        """


rule bedtools_coverage_bigwig:
    input:
        bam = "outs/{sample}/star_output/{sample}.STAR.sorted.bam",
        bamindex = "outs/{sample}/star_output/{sample}.STAR.sorted.bam.bai"
    params:
        hs1size = config["hs1size"],
    conda:
        "omics"
    resources:
        cpus_per_task = 4
    output:
        bgnotscalled = "outs/{sample}/star_output/{sample}_notscaledcov.bg",
        bg = "outs/{sample}/star_output/{sample}_cov.bg",
        bw = "outs/{sample}/star_output/{sample}_cov.bw"
    log: "logs/{sample}/deeptools_coverage_bigwig.log"
    shell:
        """
millionmapped=$(samtools idxstats {input.bam} | awk '{{sum+=$3}}END{{print sum}}') 2> {log}
mm=$(( $millionmapped / 1000000 )) 2>> {log}
scalefactor=$(bc <<<"scale=5; 1 / $mm") 2>> {log}
bedtools genomecov -ibam {input.bam} -bg > {output.bgnotscalled}.temp 2>> {log}
bedtools sort -i {output.bgnotscalled}.temp > {output.bgnotscalled} 2>> {log}
bedtools genomecov -ibam {input.bam} -bg -scale $scalefactor > {output.bg}.temp 2>> {log}
bedtools sort -i {output.bg}.temp > {output.bg} 2>> {log}
/users/mkelsey/data/tools/bedGraphToBigWig {output.bg} {params.hs1size} {output.bw} 2>> {log}
rm {output.bg}.temp 2>> {log}
rm {output.bgnotscalled}.temp 2>> {log}
        """

# rule deeptools_coverage_bedgraph:
#     input:
#         bam = "outs/{sample}/star_output/{sample}.STAR.sorted.bam",
#         bamindex = "outs/{sample}/star_output/{sample}.STAR.sorted.bam.bai"
#     conda:
#         "deeptools"
#     output:
#         coverage = "outs/{sample}/star_output/{sample}_cov.bg"
#     shell: "bamCoverage -b {input.bam} -o {output.coverage} --binSize 20 --outFileFormat bedgraph"


rule deeptools_plotcoverage:
    input:
        bam = "outs/{sample}/star_output/{sample}.STAR.sorted.bam"
    log:
        "logs/{sample}/deeptools_plotcoverage.log"
    params:
        l1hs6kb = config["l1hs6kb"],
        l1hs6kbintact = config["l1hs6kbintact"],
    conda:
        "deeptools"
    output:
        genomecoverageplot = "results/{sample}/plots/coverageGenome.pdf",
        l1hscoverageplot = "results/{sample}/plots/coverageL1HS.pdf",
        l1hsintactcoverageplot = "results/{sample}/plots/coverageL1HSintact.pdf"
    shell:
        """
plotCoverage -b {input.bam} \
--plotFile {output.genomecoverageplot} \
-n 1000000 \
--plotTitle "Whole Genome Coverage" \
--ignoreDuplicates \
--minMappingQuality 10 2> {log}

plotCoverage -b {input.bam} \
--BED {params.l1hs6kb} \
--plotFile {output.l1hscoverageplot} \
-n 1000000 \
--plotTitle "L1HS Coverage" \
--ignoreDuplicates \
--minMappingQuality 10 2>> {log}

plotCoverage -b {input.bam} \
--BED {params.l1hs6kbintact} \
--plotFile {output.l1hsintactcoverageplot} \
-n 1000000 \
--plotTitle "L1HS Coverage" \
--ignoreDuplicates \
--minMappingQuality 10 2>> {log}
        """


#note that this is coverage based and therefore not adjusted for library size
rule deeptools_plotAggcoverage:
    input:
        coverage = "outs/{sample}/star_output/{sample}_cov.bw"
    params:
        l1hs6kb = config["l1hs6kb"],
        l1hs6kbintact = config["l1hs6kbintact"]
    conda:
        "deeptools"
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



rule deeptools_plotAggcoverage2:
    input:
        coverageF = "outs/{sample}/star_output/{sample}_forward_BPM_cov.bw",
        coverageR = "outs/{sample}/star_output/{sample}_reverse_BPM_cov.bw"
    params:
        l1hs6kb = config["l1hs6kb"],
        l1hs6kbintact = config["l1hs6kbintact"]
    conda:
        "deeptools"
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

#######################################################################
#DE analysis
#######################################################################
rule mergeTElocal:
    input:
        telocal = expand("outs/{sample}/TElocal/{sample}_{{maptype}}.cntTable", sample = samples)
    params:
        sample_table = config["sample_table"]
    conda: "renv"
    log: "logs/agg/mergeTElocal_{maptype}.log"
    output:
        aggcounts = "outs/agg/TElocalCounts_{maptype}.txt",
    script:
        "scripts/mergeTEdf.R"


rule DEseq2:
    input:
        star = "outs/agg/refseq.counts.txt",
        telocal_multi = "outs/agg/TElocalCounts_multi.txt",
        telocal_uniq = "outs/agg/TElocalCounts_uniq.txt",
    params:
        sample_table = config["sample_table"],
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        levels = config["levels"],
        outputdir = "results/agg/deseq2"
    conda: "renv"
    log: "logs/agg/DEseq2.log"
    output:
        results = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"]),
        outfile = "results/agg/deseq2/outfile.txt"
    script:
        "scripts/DESeq2.R"

rule clusterprofiler:
    input:
        deresults = expand("results/agg/deseq2/star/{contrast}/results.csv", contrast = config["contrasts"]),
        normcounttable = expand("results/agg/deseq2/star/{contrast}/counttablesizenormed.csv", contrast = config["contrasts"])
    params:
        contrasts = config["contrasts"],
        genelistsforplot = config["genelistsforplot"],
        SenMayoHuman = config["SenMayoHuman"],
        sample_table = config["sample_table"],
        inputdir = "results/agg/deseq2/star",
        levels = config["levels"],
        outputdir = "results/agg/clusterprofiler"
    conda:
        "Rclusterprofiler"
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

rule gvizreduxParallel:
    input:
        sortedSTARbams = expand("outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = samples),
        DETEsbyContrast = "results/agg/repeatanalysis/allactiveDETEs.tsv"
    params:
        contrasts = config["contrasts"],
        telocaltypes = config["telocaltypes"],
        peptable = "conf/private/peptable.csv",
        refseq = config["refseq2"],
        genes = config["genes"],
        ideogram = config["ideogram"],
        intactl1s = config["intactl1s"],
        repeats = config["repeats2"],
        levels = config["levels"],
        condition_colors = config["condition_colors"],
        telocalmapping = "/users/mkelsey/data/ref/genomes/hs1/TElocal/T2T_CHM13_v2_rmsk_TE.gtf.locInd.locations",
        outputdir = "results/agg/genometracks"
    resources:
        mem_mb  = 20000,
        runtime = 30
    output:
        outfile = "temp/outfile{telocaltype}{contrast}{rtekind}{range}._{telocaltype}._{contrast}._{rtekind}._{range}.gviz.txt"
    conda:
        "repeatanalysis"
    script:
        "scripts/gviz.R"

rule repeatanalysis:
    input:
        deseq = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"]),
        telocal = expand("outs/{sample}/TElocal/{sample}_{maptype}.cntTable", sample = samples, maptype = ["multi", "uniq"])
    params:
        repeatanalysis = config["repeatanalysis"],
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        telocaltypes = config["telocaltypes"],
        levelslegendmap = config["levelslegendmap"],
        lengthreq = config["lengthreq"],
        peptable = "conf/private/peptable.csv",
        contrast_colors =config["contrast_colors"],
        condition_colors =config["condition_colors"],
        repeats = config["repeats"],
        repeatsannotatted = config["repeatsannotated"],
        telocalmapping = config["telocalmapping"],
        inputdir = "results/agg/deseq2",
        outputdir = "results/agg/repeatanalysis"
    conda:
        "repeatanalysis"
    log:
        "logs/agg/repeatanalysis.log"
    output:
        activeelementContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/activeelementContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysisactiveelementContrastplot.rst", category="repeat analysis"),
        familyContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/FamilyContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysisfamilyContrastplot.rst", category="repeat analysis"),
        combinedelementContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/CombinedContrastPlot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysiscombinedContrastplot.rst", category="repeat analysis"),
        DETEsbyContrast = "results/agg/repeatanalysis/allactiveDETEs.tsv",
        resultsdf = "results/agg/repeatanalysis/resultsdf.tsv",
        sigAluYs = "results/agg/repeatanalysis/sigAluYs.tsv",
        sigL1s = "results/agg/repeatanalysis/sigL1s.tsv",
        sigHERVKs = "results/agg/repeatanalysis/sigHERVKs.tsv",
        outfile = "results/agg/repeatanalysis/outfile.txt"
    script:
        "scripts/repeatanalysis.R"

rule derepeatReport:
    input:
        sites = "results/agg/repeatanalysis/sigL1s.tsv",
        bedgraphs = expand("outs/{sample}/star_output/{sample}_cov.bedgraph", sample = ["NT1"]),
        bigwigs = expand("outs/{sample}/star_output/{sample}_cov.bigwig", sample = ["NT1"]),
        bams = expand("outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = ["NT1", "NT2", "NT3", "PRO1", "PRO2", "PRO3"]),
        sigAluYs = "results/agg/repeatanalysis/sigAluYs.tsv",
        sigL1s = "results/agg/repeatanalysis/sigL1s.tsv",
        sigHERVKs = "results/agg/repeatanalysis/sigHERVKs.tsv"
    log: "logs/agg/derepeatReport.log"
    params:
        peptable = "conf/private/peptable.csv",
        refseq = config["refseq2"],
        genes = config["genesgz"],
        hs1fa = config["hs1sorted"],
        l1hs6kb = config["l1hs6kbintactbed"],
        ideogram = config["ideogram"],
        intactl1s = config["intactl1s"],
        repeats = config["repeats2gz"],
    conda:
        "igvreports"
    output:
        report = "results/agg/repeatanalysis/derepeatReport.html"
    shell:
        r"""
echo "[" > igvjsTrackConfig.json

# Get the array of input BAM file paths
trackdata=({input.bigwigs})

# Get the number of elements in the array
n=${{#trackdata[@]}}

# Loop over the elements of the array
for i in "${{!trackdata[@]}}"; do
  # Get the current BAM file path
  bam="${{trackdata[$i]}}"

  # Construct the JSON object for the current BAM file path
  json="{{\"autoscaleGroup\": \"A\",
         \"url\": \"$bam\"}}"

  # Append the JSON object to the output file
  if [[ $i -eq $((n-1)) ]]; then
    echo "$json" >> igvjsTrackConfig.json
  else
    echo "$json," >> igvjsTrackConfig.json
  fi
done

echo ] >> igvjsTrackConfig.json

create_report {input.sites} \
{params.hs1fa} \
--ideogram {params.ideogram} \
--sequence 1 --begin 2 --end 3 \
--flanking 10000 \
--tracks {params.genes} {params.repeats} \
--track-config igvjsTrackConfig.json \
--output {output.report} 2> {log}
        """
#\"format\": \"bam\",
         
rule repeatVariance:
    input: 
        resultsdf = "results/agg/repeatanalysis/resultsdf.tsv"
    params:
        repeatanalysis = config["repeatanalysis"],
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        telocaltypes = config["telocaltypes"],
        levelslegendmap = config["levelslegendmap"],
        lengthreq = config["lengthreq"],
        peptable = "conf/private/peptable.csv"
    output:
        corrplot1 = "results/agg/repeatanalysis/variance/{telocaltype}/{contrast}/corrplot1.pdf",
        corrplot2 = "results/agg/repeatanalysis/variance/{telocaltype}/{contrast}/corrplot2.pdf",
        corrplotcontrast = "results/agg/repeatanalysis/variance/{telocaltype}/{contrast}/corrplot{contrast}.pdf"
    notebook:
        "scripts/rteVariance.ipynb"

rule ideogram:
    input:
        DETEsbyContrast = "results/agg/repeatanalysis/allactiveDETEs.tsv"
    params:
        contrasts = config["contrasts"],
        telocaltypes = config["telocaltypes"],
        rtestoplot = config["rtestoplot"],
        sample_table = config["sample_table"],
        karyotype = config["karyotype"],
        genedensity = config["genedensity"],
        namedcolorlist = config["namedcolorlist"],
        namedmarkerlist = config["namedmarkerlist"],
        outputdir = "results/agg/repeatanalysis"
    conda:
        "repeatanalysis"
    log:
        "logs/agg/ideogram.log"
    output:
        outfile = "results/agg/repeatanalysis/ideogram_outfile.txt"
    script:
        "scripts/ideogram.R"


rule getDEelementMSA:
    input:
        "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.bed"
    params:
        consensus = "/users/mkelsey/data/ref/sequences/{rte}consensus.fa",
        hs1fa = config["hs1sorted"],
        outgroup = "/users/mkelsey/data/ref/sequences/{rte}outgroup.fa"
    conda:
        "evo"
    output:
        fa = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.fasta",
        faWithConsensus = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.withconsensus.fasta",
        faWithConsensusAndOutgroup = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.withconsensusandoutgroup.fasta",
        aln = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.aln",
        alnWithConsensus = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.WithConsensus.aln",
        alnWithConsensusAndOutgroup = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.WithConsensusAndOutgroup.aln"

    shell:
        """
if [ -s {input} ]
then
bedtools getfasta -s -fullHeader -fi {params.hs1fa} -bed {input} -fo {output.fa}
cat {params.consensus} > {output.faWithConsensus}
echo >> {output.faWithConsensus}
cat {output.fa} >> {output.faWithConsensus}
cat {params.consensus} > {output.faWithConsensusAndOutgroup} 
echo >> {output.faWithConsensusAndOutgroup}
cat {params.outgroup} >> {output.faWithConsensusAndOutgroup}
echo >> {output.faWithConsensusAndOutgroup}
cat {output.fa} >> {output.faWithConsensusAndOutgroup}
if [ $(grep ">" {output.fa} | wc -l) -gt 1 ]
then
mafft --auto {output.fa} > {output.aln}
else
touch {output.aln}
fi
mafft --auto {output.faWithConsensus} > {output.alnWithConsensus}
mafft --auto {output.faWithConsensusAndOutgroup} > {output.alnWithConsensusAndOutgroup}
else
touch {output.fa}
touch {output.faWithConsensus}
touch {output.faWithConsensusAndOutgroup}
touch {output.aln}
touch {output.alnWithConsensus}
touch {output.alnWithConsensusAndOutgroup}
fi
        """ 
#expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}.fasta", telocaltype = telocaltypes, contrast = contrasts, RTE = "L1HS", direciton = ["UP", "DOWN"])


rule evoAnalysis:
    input:
        alnWithConsensus = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}.WithConsensus.aln",
        aln = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}.aln"
    conda:
        "evo"
    params:
        basepath = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}"
    output:
        alignment = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}.alignment.pdf"
    script:
        "scripts/evoAnalysis.R"

#"results/agg/repeatanalysis/telocal_multi/condition_SEN_vs_PRO/AluY/deUP.alignment.pdf"
rule deeptools_plotAggSignal:
    input:
        coverage = expand("outs/{sample}/star_output/{sample}_cov.bw", sample = samples)
    params:
        l1hs6kb = config["l1hs6kb"],
        aluy = config["aluY"],
        hervk = config["HERVK"],
        l1hs6kbintact = config["l1hs6kbintact"],
        outputdirprofile = "results/agg/repeatanalysis",
        outputdirmatrix = "outs/agg"
    conda:
        "deeptools"
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

rule allplots:
    input:
        bamstats = "outs/{sample}/{sample}.bam.sorted.stats.txt",
        l1hsintactcoverageplot = "results/{sample}/plots/coverageL1HSintact.pdf",
        matrix = "outs/{sample}/matrix_l1hs6kb_coverage.tab.gz",
        matrixl1hs = "outs/agg/matrix_l1hs6kb_coverage.tab.gz",
        outfile3 = "results/agg/clusterprofiler/outfile.txt",
        outfile = "results/agg/repeatanalysis/genometracks/outfile.txt",
        outfile2 = "results/agg/repeatanalysis/outfile.txt",
        des2 = "results/agg/deseq2/outfile.txt"
    output:
        "allplotsDONE{sample}.txt"


#####################################################################################For the integration project only
try:
    arnasamples = peptable[peptable.batch == "alexandra"].sample_name
    marco2samples = peptable[peptable.batch == "marco2"].sample_name
    nat2019samples = peptable[peptable.batch == "nat2019"].sample_name
except:
    arnasamples = "a"
    marco2samples = "a"
    nat2019samples = "a"

### BE SURE to order samples in the input according to their order in the sample_table.csv
rule INTEGRATEfeatureCounts:
    input:
        sortedSTARbams = expand("/users/mkelsey/data/arna/outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = arnasamples) + expand("/users/mkelsey/data/senescence/outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = marco2samples) + expand("/users/mkelsey/data/marco/outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = nat2019samples)
    output:
        countsmessy = "outs/agg/INTEGRATE_refseq.counts_messy.txt",
        counts = "outs/agg/INTEGRATE_refseq.counts.txt",
        metafeaturecounts = "outs/agg/INTEGRATE_refseq.metafeature.counts.txt"
    params: 
        gtf = config['refseq']
    log: "logs/agg/featureCounts.log"
    conda:
        "deeptools"
    threads: 2
    shell: 
        """
featureCounts -p -T {threads} -t exon -a {params.gtf} -o {output.countsmessy} {input.sortedSTARbams} 2> {log}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
featureCounts -p -T {threads} -B -O -a {params.gtf} -o {output.metafeaturecounts} {input.sortedSTARbams} 2>> {log}
        """

### BE SURE to order samples in the input according to their order in the sample_table.csv
rule INTEGRATEmergeTElocal:
    input:
        telocal_multi = expand("/users/mkelsey/data/arna/outs/{sample}/TElocal/{sample}_multi.cntTable", sample = arnasamples) + expand("/users/mkelsey/data/senescence/outs/{sample}/TElocal/{sample}_multi.cntTable", sample = marco2samples) + expand("/users/mkelsey/data/marco/outs/{sample}/TElocal/{sample}_multi.cntTable", sample = nat2019samples),
        telocal_uniq =  expand("/users/mkelsey/data/arna/outs/{sample}/TElocal/{sample}_uniq.cntTable", sample = arnasamples) + expand("/users/mkelsey/data/senescence/outs/{sample}/TElocal/{sample}_uniq.cntTable", sample = marco2samples) + expand("/users/mkelsey/data/marco/outs/{sample}/TElocal/{sample}_uniq.cntTable", sample = nat2019samples)
    params:
        sample_table = config["sample_table"],
        telocaltypes = config["telocaltypes"]
    conda: "renv"
    log: "logs/agg/INTEGRATEmergeTElocal.log"
    output:
        telocal_multi = "outs/agg/INTEGRATE_TElocalCounts_multi.txt",
        telocal_uniq = "outs/agg/INTEGRATE_TElocalCounts_uniq.txt"
    script:
        "scripts/INTEGRATEmergeTEdf.R"    


rule INTEGRATEDEseq2:
    input:
        star = "outs/agg/INTEGRATE_refseq.counts.txt",
        telocal_multi = "outs/agg/INTEGRATE_TElocalCounts_multi.txt",
        telocal_uniq = "outs/agg/INTEGRATE_TElocalCounts_uniq.txt",
    params:
        peptable = "conf/private/peptable.csv",
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        levels = config["levels"],
        outputdir = "results/agg/deseq2"
    conda: "renv"
    log: "logs/agg/DEseq2.log"
    output:
        results = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"]),
        outfile = "results/agg/deseq2/INTEGRATE_outfile.txt"
    script:
        "scripts/INTEGRATEDESeq2.R"


#be sure to order the contrasts if there are multiple by the one whose naming scheme is respected in this project
rule INTEGRATEdetermineSharedDeTes:
    params:
        peptable = "conf/private/peptable.csv",
        outputdir = "results/agg/deseq2",
        telocaltypes = config["telocaltypes"],
        contraststocompare = ["condition_SEN_vs_PRO", "condition_LSEN_vs_PRO"]
    conda: "repeatanalysis"
    output:
        outfile = "outfile_sharedetes.txt"
    script:
        "scripts/INTEGRATEanalyzeDERTEs.R"


rule INTEGRATErepeatanalysis:
    input:
        deseq = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"]),
        telocal = expand("outs/agg/INTEGRATE_TElocalCounts_{maptype}.txt", sample = samples, maptype = ["multi", "uniq"])
    params:
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        telocaltypes = config["telocaltypes"],
        levelslegendmap = config["levelslegendmap"],
        peptable = "conf/private/peptable.csv",
        contrast_colors =config["contrast_colors"],
        condition_colors =config["condition_colors"],
        repeats = config["repeats"],
        telocalmapping = config["telocalmapping"],
        inputdir = "results/agg/deseq2",
        outputdir = "results/agg/repeatanalysis"
    conda:
        "repeatanalysis"
    log:
        "logs/agg/repeatanalysis.log"
    output:
        activeelementContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/activeelementContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysisactiveelementContrastplot.rst", category="repeat analysis"),
        familyContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/FamilyContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysisfamilyContrastplot.rst", category="repeat analysis"),
        combinedelementContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/CombinedContrastPlot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysiscombinedContrastplot.rst", category="repeat analysis"),
        sharedamongallcontrasts_derte = "results/agg/repeatanalysis/sharedamongallcontrasts_derte.tsv",
        DETEsbyContrast = "results/agg/repeatanalysis/allactiveDETEs.tsv",
        outfile = "results/agg/repeatanalysis/INTEGRATEoutfile.txt"
    script:
        "scripts/repeatanalysis.R"



#############################################


