log <- file(snakemake@log[[1]], open = "wt")
sink(log)
save.image()

library("readr")
library("stringr")
library("dplyr")
library("ggplot2")
library("tibble")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("eulerr")
library("ggVennDiagram")
library("Gviz")
library("GenomicRanges")
library("rtracklayer")
library("trackViewer")
library("org.Hs.eg.db")
library("httpgd")
####################################
bams <- snakemake@input[["sortedSTARbams"]]
refseq <- snakemake@params[["refseq2"]]
genes <- snakemake@params[["genes"]]
l1hs6kbintactbed <- snakemake@params[["l1hs6kbintactbed"]]
repeats <- snakemake@params[["repeats2"]]
mapper <- read.delim(snakemake@params[["telocalmapping"]], sep = "\t")
peptable <- read.csv(snakemake@params[["peptable"]])
samples <- peptable$sample_name
ideodf <- snakemake@params[["ideogram"]]
ideodf <- "/users/mkelsey/data/ref/genomes/hs1/annotations2/ideogramWithStain.bed"
levels <- c("PRO", "SEN")
levels <- snakemake@params[["levels"]]
condition_colors <- list("PRO" = "blue", "SEN" = "red")
condition_colors <- snakemake@params[["condition_colors"]]
#######################################
# getOption("Gviz.scheme")
# ## [1] "default"
scheme <- getScheme()
scheme$GeneRegionTrack$arrowHeadMaxWidth <- 5
scheme$GeneRegionTrack$arrowHeadWidth <- 5
scheme$GeneRegionTrack$col <- "black"
scheme$GeneRegionTrack$fill <- "#FFD58A"
addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")

# scheme$GeneRegionTrack$shape <- c("smallArrow", "arrow")


tab <- table(mapper$TE)
getPos <- function(rte, tab) {
    if (unname(tab[rte]) > 1) {
        return(NULL)
    } else {
        rterow <- mapper %>% filter(TE == rte)
        val <- rterow["chromosome.start.stop.strand"]
        chr <- str_split(val, ":")[[1]][1]
        coords <- str_split(val, ":")[[1]][2]
        strand <- str_split(val, ":")[[1]][3]
        start <- str_split(coords, "-")[[1]][1]
        stop <- str_split(coords, "-")[[1]][2]

        return(list(chr = chr, start = as.numeric(start), stop = as.numeric(stop), strand = strand))
    }
}

genome <- "hs1"
ideodf <- read.delim(ideodf, header = TRUE)
strack <- SequenceTrack("/users/mkelsey/data/ref/genomes/hs1/hs1.sorted.fa")
# localrtetrack <- AnnotationTrack(start = start, width = stop - start, chromosome = chr, strand = strand, id = rte, genome = genome)
ideogramtrack <- IdeogramTrack(genome = genome, bands = ideodf)

# fix! levels <- samples


bamsubset <- bams[(grepl("PRO", bams) | grepl("NT", bams))]

svec <- bamsubset
alignmentTrackList <- list()
for (i in seq(length(svec))) {
    sample <- str_split_1(basename(svec[i]), "\\.")[1]
    condition <- peptable[peptable$sample_name == sample, ]$condition
    fillcolor <- condition_colors[[condition]]
    assign(
        paste0("alTrack", i),
        AlignmentsTrack(svec[i],
            name = sample,
            isPaired = TRUE,
            genome = "hs1",
            ylim = c(0, 25),
            type = "coverage",
            fill = fillcolor
        )
    )
    track <- get(paste0("alTrack", i))
    alignmentTrackList[sample] <- track
}

# genesdf = read.delim(genes, header = FALSE, sep = "")
# names(genesdf)[1] = "chromosome"
# names(genesdf)[3] = "feature"
# names(genesdf)[4] = "start"
# names(genesdf)[5] = "end"
# names(genesdf)[7] = "strand"
# names(genesdf)[16] = "gene"

genesTrack <- GeneRegionTrack(genes,
    genome = genome,
    showFeatureId = TRUE,
    showID = TRUE,
    name = "RefSeqCurated",
    stacking = "squish",
    collapseTranscripts = TRUE,
    shape = "arrow"
)

repeatTrack <- GeneRegionTrack(repeats,
    name = "RepeatMasker",
    genome = genome,
    showFeatureId = TRUE,
    showID = TRUE,
    stacking = "squish",
    shape = "arrow",
    fill = "#95c8ff"
)

dertelist <- c("L1HS_dup291")
for (rte in dertelist) {
    location <- getPos(rte, tab)
    genome <- "hs1"
    chr <- location["chr"][[1]]
    strand <- location["strand"][[1]]
    start <- location["start"][[1]]
    stop <- location["stop"][[1]]

    length <- stop - start
    gtrack <- GenomeAxisTrack(name = paste0(length, "bp"))

    tracks <- c(gtrack, alignmentTrackList, genesTrack, repeatTrack)
    lengthtracks <- length(tracks)
    ht <- HighlightTrack(
        trackList = tracks,
        start = start, end = stop,
        chromosome = chr,
        inBackground = TRUE,
        fill = "#d5e2ff",
        col = "transparent"
    )

    pdf("z.pdf", width = 6, height = 8)
    pl <- plotTracks(
        c(ideogramtrack, ht),
        chromosome = chr,
        from = start,
        to = stop,
        extend.right = 5000,
        extend.left = 5000,
        sizes = c(0.2, rep(0.3, lengthtracks)),
        showFeatureId = TRUE,
        featureAnnotation = "id",
        shape = "fixedArrow",
        collapseTranscripts = TRUE,
        showId = TRUE,
        col = "black",
        fontcolor = "black",
        fontcolor.group = "black",
        background.panel = NULL,
        background.title = "#dadada",
        col.title = "black",
        col.border.title = "black",
        col.axis = "black",
        col = "black",
        showSampleNames = TRUE
    )
    dev.off()
}

#####################

x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)

dev.off()
