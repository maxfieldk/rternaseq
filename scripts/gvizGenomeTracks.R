log <- file(snakemake@log[[1]], open="wt")
sink(log)

library('readr')
library("stringr")
library('dplyr')
library("ggplot2")
library('tibble')
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("eulerr")
library("ggVennDiagram")
library('Gviz')
library(GenomicRanges)
library('rtracklayer')
library('trackViewer')
library("org.Hs.eg.db")
####################################
save.image()

bams = snakemake@input[["sortedbamSTAR"]]
refseq = snakemake@params[["refseq"]]
l1hs6kbintactbed = snakemake@params[["l1hs6kbintactbed"]]
repeatsbed = snakemake@params[["repeatsbed"]]
mapper = read.delim(snakemake@params[["telocalmapping"]], sep = "\t")
#######################################

tab = table(mapper$TE)
getPos <- function(rte, tab) {
    if (unname(tab[rte]) > 1) {
        return(NULL)
    } else {
        rterow = mapper %>% filter(TE == rte)
        val = rterow["chromosome.start.stop"]
        chr = str_split(val, ":")[[1]][1]
        coords = str_split(val, ":")[[1]][2]
        strand = str_split(val, ":")[[1]][3]
        start = str_split(coords, "-")[[1]][1]
        stop = str_split(coords, "-")[[1]][2]
        
        return(list(chr = chr, start = as.numeric(start),stop = as.numeric(stop), strand = strand))
    }
}

gtrack <- GenomeAxisTrack()
strack = SequenceTrack("/users/mkelsey/data/ref/genomes/hs1/hs1.sorted.fa")
repeattrack = AnnotationTrack(repeatsbed, name = "repeat masker")


levels = snakemake@params["samples"]

alignmenttracklist = list()
datatracklist = list()
for (i in seq(length(bams))) {
    assign(paste0("alTrack", i), AlignmentsTrack(bams[i],name = basename(bams[i]), isPaired = TRUE))
    track = get(paste0("alTrack", i))
    alignmenttracklist[basename(bams[i])] = track

    assign(paste0("dataTrack", i),  DataTrack(range = bams[i], genome = "hs1", type = "l", 
                name = "coverage",  groups = factor(levels[[i]], levels = levels), legend = TRUE))
    track = get(paste0("dataTrack", i))
    datatracklist[basename(bams[i])] = track

}

dertelist = c("L1HS_dup28")
for (rte in dertelist) {
    location = getPos(rte, tab)
    
    localrtetrack <- AnnotationTrack(start = location["start"][[1]], width = location["stop"][[1]] - location["start"][[1]], chromosome = location["chr"][[1]], strand = location['strand'][[1]], id = rte)
    pdf("scripts/zzz.pdf", width = 6, height = 18)
    plotTracks(c(strack, gtrack, localrtetrack, alignmenttracklist),
    chromosome = coords["chr"][[1]],
    from = coords["start"][[1]],
    to = coords["stop"][[1]], extend.left = 0.1, extend.right = 0.1,  ylim = c(0, 50),
    featureAnnotation = "id")

    reps = 3

    length(datatracklist)
    counter = 1
    for (i in seq(length(datatracklist))) {
        if 
    }
    ot <- OverlayTrack(trackList=list(dTrack2, dTrack3))
    ylims <- extendrange(range(c(values(dTrack3), values(dTrack2))))


    lim = c(1001000, 1001500)
    pdf("scripts/zz.pdf", width = 6, height = 6)
    plotTracks(list(gtrack, localrtetrack, ot),chromosome = "chr19", from = lim[1], to = lim[2], ylim = c(0, 4000))
    dev.off()

dev.off()


}



















pdf("scripts/z.pdf", width = 6, height = 20)
plotTracks(list(strack,gtrack, alTrack1, alTrack2), chromosome = "chr19",
    from= 1001450, to = 1001550)


dev.off()




#####################

x <- data.frame()
write.table(x, file=snakemake@output[['outfile']], col.names=FALSE)


