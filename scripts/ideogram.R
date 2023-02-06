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
library("RIdeogram")
library("httpgd")


####### functions
get_ideogram_labels <- function(
    rtedf,
    elementlist,
    directions = c("UP", "DOWN"),
    contrast = "condition_LSEN_vs_PRO",
    tecountstype = "telocal_multi",
    namedmarkerlist = c("L1" = "circle", "Alu" = "box", "ERVK" = "triange", "L1HS" = "circle", "AluY" = "box", "HERVK-int" = "triangle"),
    namedcolorlist = c(
        "L1UP" = "ff0000", "AluUP" = "ff9d00", "ERVKUP" = "ff7be2", "L1HSUP" = "ff0000", "AluYUP" = "ff9d00", "HERVK-intUP" = "ff7be2",
        "L1DOWN" = "3846ff", "AluDOWN" = "0be7e4", "ERVKDOWN" = "8c00ff", "L1HSDOWN" = "3846ff", "AluYDOWN" = "0be7e4", "HERVK-intDOWN" = "8c00ff"
    )) {
    tempdf <- rtedf[, c(1, 3, 4, 5, 7)]
    colnames(tempdf) <- c("Kind", "Chr", "Start", "End", "Direction")
    print("a")
    print(head(tempdf))
    tempdf <- tempdf %>% filter(Kind %in% elementlist)
    print("a1")
    print(head(tempdf))
    tempdf <- tempdf %>%
        mutate(Shape = namedmarkerlist[Kind]) %>%
        mutate(color = namedcolorlist[paste0(Kind, Direction)]) %>%
        mutate(Type = paste(Kind, Direction, sep = " "))
    # now the df is in a format that RIdeogram can read
    df <- tempdf[, c("Type", "Shape", "Chr", "Start", "End", "color")]
    print("b")
    print(head(df))
    return(df)
}

makeIdeogram <- function(
    #' Be sure to have svgsavepath end kn .svg
    rtedf,
    elementlist,
    karyotype,
    gff,
    directions = c("UP", "DOWN"),
    svgsavepath = "chromosome.svg",
    contrast = "condition_LSEN_vs_PRO", tecountstype = "telocal_multi",
    namedmarkerlist = c("L1" = "circle", "Alu" = "box", "ERVK" = "triange", "L1HS" = "circle", "AluY" = "box", "HERVK-int" = "triangle"),
    namedcolorlist = c("L1" = "ff0000", "Alu" = "ffa500", "ERVK" = "0000ff", "L1HS" = "ff0000", "AluY" = "ffa500", "HERVK-int" = "0000ff")) {
    karyotypedf <- read.delim(karyotype, sep = "\t")
    gene_density <- GFFex(
        input = gff,
        karyotype = karyotype,
        feature = "exon",
        window = 1000000
    )
    dertes <- read.delim(rtedf, sep = "\t", header = FALSE)
    rownames(dertes) <- 1:nrow(dertes)

    labels <- get_ideogram_labels(dertes, elementlist,
        directions = directions,
        contrast = contrast,
        tecountstype = tecountstype
    )
    ideogram(
        karyotype = karyotypedf,
        overlaid = gene_density,
        label = labels,
        label_type = "marker",
        output = svgsavepath
    )
    convertSVG(svgsavepath, gsub(".svg", ".png", svgsavepath), device = "png")
}

########
rtedf <- "/users/mkelsey/data/marco/results/agg/repeatanalysis/allactiveDETEs.tsv"
gff <- "/users/mkelsey/data/ref/genomes/hs1/annotations/hs1.110.20220412.ncbiRefSeqUCSCstyle.gtf"
karyotype <- "/users/mkelsey/data/ref/genomes/hs1/karyotype.tsv"
rtenames <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int")
elementlist <- rtenames[4]
savepath <- "chromosometwo.svg"

elementlist
directions <- c("UP", "DOWN")



makeIdeogram(rtedf, c("L1HS"), karyotype, gff, svgsavepath = "ZZZZZ.svg")






directions <- c("UP", "DOWN")
contrast <- "condition_LSEN_vs_PRO"
tecountstype <- "telocal_multi"
namedmarkerlist <- c("L1" = "circle", "Alu" = "box", "ERVK" = "triange", "L1HS" = "circle", "AluY" = "box", "HERVK-int" = "triangle")
namedcolorlist <- c(
    "L1UP" = "#8c00ff", "AluUP" = "#ff9d00", "ERVKUP" = "#ff7be2", "L1HSUP" = "8c00ff", "AluYUP" = "ff9d00", "HERVK-intUP" = "ff7be2",
    "L1DOWN" = "#00ff62", "AluDOWN" = "#0be7e4", "ERVKDOWN" = "#8c00ff", "L1HSDOWN" = "00ff62", "AluYDOWN" = "0be7e4", "HERVK-intDOWN" = "8c00ff"
)

dertes <- read.delim(rtedf, sep = "\t", header = FALSE)
tempdf <- dertes[, c(1, 3, 4, 5, 7)]
colnames(tempdf) <- c("Type", "Chr", "Start", "End", "Direction")
tempdf <- tempdf %>%
    filter(Type %in% elementlist) %>%
    filter(Direction %in% directions)
tempdf <- tempdf %>%
    mutate(Shape = namedmarkerlist[Type]) %>%
    mutate(color = namedcolorlist[paste0(Type, Direction)])
# now the df is in a format that RIdeogram can read
df <- tempdf[, c("Type", "Shape", "Chr", "Start", "End", "color")]


karyotypedf <- read.delim(karyotype, sep = "\t")
gene_density <- GFFex(
    input = gff,
    karyotype = karyotype,
    feature = "exon",
    window = 1000000
)


svgsavepath <- "testtest.svg"
labels <- df
ideogram(
    karyotype = karyotypedf,
    overlaid = gene_density,
    label = labels,
    label_type = "marker",
    output = svgsavepath
)
convertSVG(svgsavepath, gsub(".svg", ".png", svgsavepath), device = "png")
