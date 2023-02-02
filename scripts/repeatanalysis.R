#note! there are dataset specific parameters set! need to edit accordingly

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
save.image()

#order matters for the colors!
contrast_colors = snakemake@params[["contrast_colors"]]
condition_colors = snakemake@params[["condition_colors"]]
######## MAIN FUNCTIONS
####PLOTTING
plotRTE = function(rte, df, lims= c(-1, 11), title = "", colors = c("red", "blue", "grey"), alpha = 1, number_to_sample = NULL) {
    matches = grep("ERV.:LTR", results$Geneid)
    dat = results[matches,]
    if (is.null(number_to_sample)) {
        rteplot = dat %>% ggplot() +
            geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance)) +
            geom_abline(intercept = 0)+ 
            coord_fixed() +
            xlim(lims) + ylim(lims) + 
            scale_color_manual(values = colors) +
            theme_cowplot() +
            ggtitle(title) +
            panel_border(color = "black", linetype = 1, remove = FALSE) +
            theme(axis.line=element_blank()) +
            labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))
        return(rteplot)
    } else {
        number_to_sample = min(length(rownames(dat)), number_to_sample)
        datasampled = dat[sample(rownames(dat), number_to_sample,replace=FALSE),]
        rteplot = datasampled %>% ggplot() +
            geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance)) +
            geom_abline(intercept = 0)+ 
            coord_fixed() +
            xlim(lims) + ylim(lims) + 
            scale_color_manual(values = colors) +
            theme_cowplot() +
            ggtitle(title) +
            panel_border(color = "black", linetype = 1, remove = FALSE) +
            theme(axis.line=element_blank()) +
            labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))
        return(rteplot)

    }
    
}

plotAggRTE = function(df, coltogroupby = "ActiveTE", valuescol = "log2FoldChange") {
    agg = df %>% filter(.data[[coltogroupby]] != "Other") %>% ggplot(aes(x = .data[[coltogroupby]], y = .data[[valuescol]])) +
        geom_violin(draw_quantiles = c(0.5)) +
        stat_summary(fun = "mean",geom = "point",color = "black") +
        geom_hline(aes(yintercept = 0), color = 'red') +
        ggtitle("RTEs") +
        coord_fixed() +
        ylim(c(-11,11)) +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        theme(aspect.ratio=1)
    return(agg)
}

euPlot = function(fit, fill = c("red", "blue"), legend = TRUE, labels = FALSE, main = "") {
    plot(
        fit,
        main = main,
        fill = fill,
        legend = legend,
        quantities = TRUE,
        labels = labels)
}
######
##########################
contrasts = snakemake@params[["contrasts"]]
levelslegendmap =snakemake@params[["levelslegendmap"]]

for (counttype in snakemake@params[["telocaltypes"]]) {

for (contrast in contrasts) {
    #######DATA WRANGLING
    contrast_level_2 = unlist(strsplit(contrast, "_", fixed = TRUE))[[2]]
    contrast_base_level = unlist(strsplit(contrast, "_", fixed = TRUE))[[4]]
    conditions = c(contrast_base_level, contrast_level_2)

    ddsres = read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast, "results.csv", sep = "/"))
    ddscounts = read_csv(paste(snakemake@params[["inputdir"]], counttype,  contrast, "counttablesizenormed.csv", sep = "/"))
    ddsrlogcounts = read_csv(paste(snakemake@params[["inputdir"]], counttype,  contrast, "rlogcounts.csv", sep = "/"))
    
    outputdir = paste(snakemake@params[["outputdir"]], counttype, sep = "/")

    colnames(ddsres)[1] <- "Geneid"
    colnames(ddscounts)[1] <- "Geneid"
    colnames(ddsrlogcounts)[1] <- "Geneid"

    sample_table = read.csv(snakemake@params[["sample_table"]])
    
    avoidzero = 1
    meancols = c()
    log2meancols = c()
    rlogmeancols = c()

    for (condition in conditions) {
        condition_samples = filter(sample_table, condition == {{condition}})$sample_name
        s = ""
        for (e in condition_samples) {
            if (s == "") {
                s = paste0("`", e, "`")
            } else { 
                s = paste0(s, "+", "`", e, "`")
            }
        }
        s = paste0("(", s, ")", "/", length(condition_samples))
        meancol = paste0(condition, "mean")
        ddscounts = ddscounts %>% mutate({{meancol}} := eval(parse(text=s)))
        rlogmeancol = paste0("rlog", condition, "mean")
        ddsrlogcounts = ddsrlogcounts %>% mutate({{rlogmeancol}} := eval(parse(text=s)))
        log2meancol = paste0("log2",condition, "mean")
        ddscounts = ddscounts %>% mutate({{log2meancol}} := log2(.data[[meancol]]) )
        meancols = c(meancols, meancol)
        log2meancols = c(log2meancols, log2meancol)
        rlogmeancols = c(rlogmeancols, rlogmeancol)
    }

    columns_to_retain = c("Geneid",meancols, log2meancols)
    rlogcolumns_to_retain = c("Geneid", rlogmeancols)

    log2sub = ddscounts[,columns_to_retain]
    rlogsub = ddsrlogcounts[,rlogcolumns_to_retain]

    dflist = list(ddsres,log2sub, rlogsub)

    results = Reduce(function(x, y) merge(x, y, by="Geneid"), dflist, accumulate=FALSE)
    results = results %>% mutate(Significance = ifelse(padj < 0.05, ifelse(padj < 0.001, "< 0.001", "< 0.05"), "> 0.05"))
    genes = results$Geneid

    results["Family"] = "Other"
    matches = grep("L1:LINE", genes)
    results[matches,"Family"] = "L1"
    matches = grep("Alu:SINE", genes)
    results[matches,"Family"] = "Alu"
    matches = grep("ERVK:LTR", genes)
    results[matches,"Family"] = "HERV"

    results["ActiveTE"] = "Other"
    matches = grep("L1HS:L1", genes)
    results[matches,"ActiveTE"] = "L1HS"
    matches = grep("AluY:Alu", genes)
    results[matches,"ActiveTE"] = "AluY"
    matches = grep("HERVK", genes)
    results[matches,"ActiveTE"] = "HERVK"

    ################################################PLOTING
    ##### plot settings!

    #quanttype = "log2"
    quanttype = "rlog"
    xval = paste0(quanttype, contrast_base_level, "mean")
    yval = paste0(quanttype, contrast_level_2, "mean")
    l1_lims = c(-1, 11)
    alu_lims = c(-1, 11)
    herv_lims = c(-1,11)

    ##### ACTIVE RTE PLOTS
    l1hs = plotRTE("L1HS:L1", results, l1_lims, title = "L1HS",)
    aluy = plotRTE("AluY:Alu", results, alu_lims, title = "AluY",)
    hervk = plotRTE("HERVK(.)*int", results, herv_lims, title = "HERVK",)
    aggactive = plotAggRTE(results, coltogroupby = "ActiveTE")

    legend <- get_legend(l1hs + theme(legend.box.margin = margin(0, 0, 0, 1)))
    p = plot_grid(aggactive + theme(legend.position="none"),
            l1hs + theme(legend.position="none"),
            aluy + theme(legend.position="none"),
            hervk + theme(legend.position="none"), nrow = 1,
            rel_widths = c(1,1,1,1), labels = "AUTO",
            align = "vh",
            axis = 'bt')

    pdf(paste(outputdir, contrast, paste0("activeelementContrastplot", ".pdf"), sep = '/'), width=14, height=4)
    print(plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep=' ')))
    dev.off()

    ##### FAMILY LEVEL RTE PLOTS
    ##settings for all plots
    number_to_sample = 1000
    alpha = 0.9

    l1 = plotRTE("L1:LINE", results, l1_lims, title = "L1", alpha = alpha, number_to_sample = number_to_sample)
    alu = plotRTE("Alu:SINE", results, alu_lims,title = "Alu", alpha = alpha, number_to_sample = number_to_sample)
    herv = plotRTE("ERV.:LTR", results, herv_lims,title = "ERV", alpha = alpha, number_to_sample = number_to_sample)
    aggfamily = plotAggRTE(results, coltogroupby = "Family")
    legend <- get_legend(l1 + theme(legend.box.margin = margin(0, 0, 0, 1)))
    p = plot_grid(aggfamily + theme(legend.position="none"),l1 + theme(legend.position="none"), alu +
            theme(legend.position="none"), herv +
            theme(legend.position="none"), nrow = 1,
            rel_widths = c(1,1,1,1), labels = "AUTO",
            align = "vh",
            axis = 'bt')
    pdf(paste(outputdir, contrast, paste0("familyContrastplot", ".pdf"), sep = '/'), width=14, height=4)
    print(plot_grid(p, legend, rel_widths = c(3, .4))+ ggtitle(paste("DE repeats", contrast, sep=' ')))
    dev.off()

    ########## Plot activete and family together
    legend <- get_legend(l1 + theme(legend.box.margin = margin(0, 0, 0, 1)))

    p = plot_grid(aggfamily + theme(legend.position="none"),
            l1 + theme(legend.position="none"),
            alu + theme(legend.position="none"),
            herv + theme(legend.position="none"),
            aggactive + theme(legend.position="none"),
            l1hs + theme(legend.position="none"), 
            aluy + theme(legend.position="none"),
            hervk + theme(legend.position="none"),
            nrow = 2,
            rel_widths = c(1,1,1,1), 
            labels = "AUTO",
            align = "vh",
            axis = 'bt')
    pdf(paste(outputdir, contrast, paste0("combinedContrastplot", ".pdf"), sep = '/'), width=14, height=8)
    print(plot_grid(p, legend, rel_widths = c(3, .4))+ ggtitle(paste("DE repeats", contrast, sep=' ')))
    dev.off()
################################################END#PLOTING
}
}

# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:

# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. Customization of this outlier filtering and description of functionality for replacement of outlier counts and refitting is described below
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. Description and customization of independent filtering is described below


################################################Venn diagrams

for (counttype in snakemake@params[["counttypes"]]) {
    outputdir = paste(snakemake@params[["outputdir"]], counttype, sep = "/")
    UP = list()
    DOWN = list()

    eulerr_options(quantities = list(font = 2, cex = 1.5), legend = list(side = "bottom", cex = 1.5))

    for (contrast in contrasts) {
        df = read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast , "results.csv", sep = "/"))
        colnames(df)[1] <- "Geneid"
        dfUP = df[df$pvalue < 0.05 & df$log2FoldChange > 0, 'Geneid'] %>% na.omit()
        dfDOWN = df[df$pvalue < 0.05 & df$log2FoldChange < 0, 'Geneid'] %>% na.omit()
        UP[[gsub("condition_", "", contrast)]] = dfUP
        DOWN[[gsub("condition_", "", contrast)]] = dfDOWN
    }
    fit <- euler(UP, shape = "ellipse")
    upplot <- euPlot(fit)

    fit <- euler(DOWN, shape = "ellipse")
    downplot <- euPlot(fit)
    
    plots = list(upplot, downplot)
    grid = plot_grid(plotlist = plots, labels = c("UP", "DOWN"), ncol = 1)
    pdf(paste(outputdir, paste0("VennDiagram", "allDEGS", ".pdf"), sep = '/'), width=4, height=6)
    print(grid)
    dev.off()
}
####################
master_deRTEs = list()
for (counttype in snakemake@params[["telocaltypes"]]) {
    outputdir = paste(snakemake@params[["outputdir"]], counttype, sep = "/")
    pattern = list("L1:LINE", "Alu:SINE", "ERVK:LTR", "L1HS:L1", "AluY:Alu", "HERVK(.)*int")
    names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int")
    plots = list()
    deRTEs = list()
    for (e in names(pattern)) {
        rte = e
        contrastL = list()
        for (contrast in contrasts) {
            df = read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast , "results.csv", sep = "/"))
            colnames(df)[1] <- "Geneid"
            matches = grep(pattern[e], df$Geneid)
            searched_element = df[matches,]
            sigsearched = list(searched_element[searched_element$padj < 0.05, 'Geneid'] %>% na.omit())  
            names(sigsearched) <- gsub("condition_", "", contrast)
            contrastL = c(contrastL, sigsearched)
            }
        deRTEs[[e]] = contrastL
        fit <- euler(contrastL, shape = "ellipse")
        p = euPlot(fit, contrast_colors, main = rte)
        pdf(paste(outputdir, paste0("VennDiagram", rte, ".pdf"), sep = '/'), width=4, height=4)
        print(p)
        dev.off()

        p_no_title = euPlot(fit, contrast_colors, legend = TRUE)
        plots = c(plots, list(p_no_title))
    }
master_deRTEs[[counttype]] = deRTEs
grid = plot_grid(plotlist = plots, labels = names(pattern), ncol = 3)
pdf(paste(outputdir, paste0("VennDiagram", "RTEs", ".pdf"), sep = '/'), width=9, height=6)
print(grid)+ ggtitle(paste("DE repeats", contrast, sep=' '))
dev.off()
}

master_deRTEs_df <- as.data.frame(do.call(cbind, master_deRTEs))
try(write.table(master_deRTEs_df, file=snakemake@output[['deRTEsDF']]), TRUE)

########### Are de RTEs the same?

#this will be: for each counttype, what RTEs were DE in all contrasts 
des_by_counttype = list()
for (counttype in snakemake@params[["telocaltypes"]]) {
    outputdir = paste(snakemake@params[["outputdir"]], counttype, sep = "/")
    pattern = list("L1:LINE", "Alu:SINE", "ERVK:LTR", "L1HS:L1", "AluY:Alu", "HERVK(.)*int")
    names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int")
    des = list()
    for (e in names(pattern)) {
        decontrast = list()
        for (contrast in contrasts) {
            decontrast = c(decontrast, list(master_deRTEs[[counttype]][[e]][[gsub("condition_", "", contrast)]]))
        }
    set_des = Reduce(intersect, decontrast)
    des[[e]] = set_des
    }
    des_by_counttype[[counttype]]  = des
}

#this will be: for each RTE, which were DE in all contrasts AND for all count methods  
des_by_rte = list()
for (e in names(pattern)) {
    des = list()
    for (countype in snakemake@params[["telocaltypes"]]) {
        des = c(des, list(des_by_counttype[[counttype]][[e]]))
    }
    set_des = Reduce(intersect, des)
    des_by_rte[[e]] = set_des
}



# ############################### GVIZ
# outputdir = paste(snakemake@params[["outputdir"]], "genometracks", sep = "/")

# bams = snakemake@input[["sortedbamSTAR"]]
# refseq = snakemake@params[["refseq"]]
# l1hs6kbintactbed = snakemake@params[["l1hs6kbintactbed"]]
# repeatsbed = snakemake@params[["repeatsbed"]]
# hs1sorted = snakemake@params[["hs1sorted"]]
# mapper = read.delim(snakemake@params[["telocalmapping"]], sep = "\t")

# #### functions
# getPos <- function(rte, tab) {
#     if (unname(tab[rte]) > 1) {
#         return(NULL)
#     } else {
#         rterow = mapper %>% filter(TE == rte)
#         val = rterow["chromosome.start.stop.strand"]
#         chr = str_split(val, ":")[[1]][1]
#         coords = str_split(val, ":")[[1]][2]
#         strand = str_split(val, ":")[[1]][3]
#         start = str_split(coords, "-")[[1]][1]
#         stop = str_split(coords, "-")[[1]][2]
        
#         return(list(chr = chr, start = as.numeric(start),stop = as.numeric(stop), strand = strand))
#     }
# }

# #####
# tab = table(mapper$TE)
# gtrack = GenomeAxisTrack()
# strack = SequenceTrack(hs1sorted)
# repeattrack = AnnotationTrack(repeatsbed, name = "repeat masker")




# levels = snakemake@params["samples"]

# alignmenttracklist = list()
# datatracklist = list()
# for (i in seq(length(bams))) {
#     assign(paste0("alTrack", i), AlignmentsTrack(bams[i],name = basename(bams[i]), isPaired = TRUE))
#     track = get(paste0("alTrack", i))
#     alignmenttracklist[basename(bams[i])] = track

#     assign(paste0("dataTrack", i),  DataTrack(range = bams[i], genome = "hs1", type = "l", name = basename(bams[i]),  groups = factor(levels[i], levels = levels), legend = TRUE))
#     track = get(paste0("dataTrack", i))
#     datatracklist[basename(bams[i])] = track

# }
# for (name in names(des_by_rte)) {
#     dertelist = des_by_rte[name]
#     for (rte in unlist(dertelist)) {
#         rte = str_split(rte, ":")[[1]][1]
#         print(rte)
#         location = getPos(rte, tab)
#         chr = location["chr"][[1]]
#         start = location["start"][[1]]
#         stop = location["stop"][[1]]
#         strand = location["strand"][[1]]
        
#         localrtetrack <- AnnotationTrack(start = start, width = stop - start, chromosome = chr, strand = strand, id = rte)
#         dir = paste(outputdir, "consitentlyDE", sep = '/')
#         if (!dir.exists(dir)){
#             dir.create(dir)}
#         pdf(paste(outputdir, "consitentlyDE", paste0(name, "_", rte, ".pdf"), sep = '/'), width=9, height=50)
#         plotTracks(c(strack, gtrack, localrtetrack, alignmenttracklist),
#         chromosome = chr,
#         from = start,
#         to = stop,
#         extend.left = 0.1, extend.right = 0.1,  ylim = c(0, 50),
#         featureAnnotation = "id")
#         dev.off()

# dt1 = DataTrack(range = "/users/mkelsey/data/marco/outs/SRR6515349/SRR6515349_cov.bw", 
#         chromosome = chr,
#         from = start,
#         to = stop, type = "histogram", name = basename(bams[1]), legend = TRUE)

# dt2 = AlignmentsTrack(range = bams[2],
#         chromosome = chr,
#         from = start,
#         to = stop,genome = "hs1", type = "coverage", name = basename(bams[2]), legend = TRUE)



#  ylims <- extendrange(range(c(values(dt1), values(dt2))))
#  range(c(values(dTrack3), values(dTrack2))))


#         pdf("cooltime.pdf", width=9, height=4)
#         plotTracks(c(strack, gtrack, dt1, dt2),
#         chromosome = chr,
#         from = start,
#         to = stop,
#         extend.left = 0.1, extend.right = 0.1)
#         dev.off()
        
#         reps = 3
#         if ((length(datatracklist) %% reps) == 0) {
#             timestocycle = length(datatracklist) / reps

#             tracks_grouped = list(strack, localrtetrack)
#             for (i in seq(timestocycle)) {
#                 i = i-1
#                 tracks = datatracklist[((reps*i)+1):((reps*i)+3)]
#                 ot <- OverlayTrack(trackList=tracks)
#                 tracks_grouped = c(tracks_grouped, ot)
#             }
        
#             pdf(paste(outputdir, "consitentlyDE", paste0(name, "_", rte, "Grouped", ".pdf"), sep = '/'), width=9, height=6)
#             plotTracks(tracks_grouped,
#                 chromosome = chr,
#                 from = start,
#                 to = stop,
#                 extend.left = 0.1, extend.right = 0.1,  ylim = c(0, 1000),
#                 featureAnnotation = "id")
#             dev.off()
#         }
        
#     }

# }


# getOption("Gviz.scheme")
# ## [1] "default"
# scheme <- getScheme()
# scheme$background.title = "darkblue"
# scheme$background.panel = "#FFFEDB"
# addScheme(scheme, "myScheme")
# options(Gviz.scheme = "myScheme")



# tracks_grouped = list(strack, localrtetrack)
# for (i in seq(timestocycle)) {
#     i = i-1
#     tracks = datatracklist[((reps*i)+1):((reps*i)+3)]
#     ot <- OverlayTrack(trackList=tracks)
#     tracks_grouped = c(tracks_grouped, ot)
# }

# datatracklist


# pdf(paste("play.pdf", sep = '/'), width=9, height=6)
# plotTracks(tracks_grouped[4],
#     chromosome = chr,
#     from = start,
#     to = stop,
#     extend.left = 0.1, extend.right = 0.1,  ylim = c(0, 1000),
#     featureAnnotation = "id")
# dev.off()


# pdf(paste("play.pdf", sep = '/'), width=9, height=6)
# plotTracks(ot,
#     chromosome = chr,
#     from = start,
#     to = stop,
#     type = c("smooth", "p"),
#     each = 3,
#     extend.left = 0.1, extend.right = 0.1,  ylim = c(0, 1000),
#     featureAnnotation = "id")
# dev.off()














#################### RTE in genome histograms
outputdir = snakemake@params[["outputdir"]]
repeats = snakemake@params[["repeats"]]
repmasker = read_delim(repeats, col_names=FALSE) %>% mutate_at(c('X2', 'X3'), as.numeric) %>% mutate(length = X3-X2)


pattern = c("LINE/L1", "SINE/Alu", "LTR/ERVK", "L1HS", "AluY", "HERVK(.)*int")
names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int")
activeLengthReq = c(1,0,0,6000,290,7000)
maxlength = c(7000, 400, 15000,7000, 400, 15000 )
colors = c("deeppink4", "deepskyblue4", "darkorange4", "deeppink1", "deepskyblue1", "darkorange1")

plots = list()

for (i in seq(length(pattern))) {
print(i)
print(pattern[i])
p = repmasker %>% filter(length < maxlength[i]) %>%  filter( grepl(pattern[i], X4)) %>% ggplot() +
geom_histogram(aes(x=length), fill = colors[i]) + ggtitle(names(pattern)[i])+ theme_cowplot()
plots = c(plots, list(p))
}

grid = plot_grid(plotlist = plots, labels = "AUTO")
pdf(paste(outputdir, paste0("histogram", "allRTEs", ".pdf"), sep = '/'), width=10, height=6)
print(grid)
dev.off()

##############

x <- data.frame()
write.table(x, file=snakemake@output[['outfile']], col.names=FALSE)

#############


 