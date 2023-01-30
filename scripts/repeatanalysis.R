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
library('rtracklayer')
library('trackViewer')
library("org.Hs.eg.db")
save.image()
levelslegendmap =snakemake@params[["levelslegendmap"]]

for (counttype in snakemake@params[["telocaltypes"]]) {

for (contrast in snakemake@params[["contrasts"]]) {
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

    #ddsrlogcounts = ddsrlogcounts %>% mutate(rlogprolmean = ((SRR6515349 +SRR6515350 +SRR6515351)/3), rlogesenmean = ((SRR6515352 +SRR6515353 +SRR6515354)/3), rloglsenmean = ((SRR6515355 +SRR6515356 +SRR6515357)/3) )
    #ddscounts = ddscounts %>% mutate(prolmean = ((SRR6515349 +SRR6515350 +SRR6515351)/3), esenmean = ((SRR6515352 +SRR6515353 +SRR6515354)/3), lsenmean = ((SRR6515355 +SRR6515356 +SRR6515357)/3) )
    
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

    #plot settings!
    #quanttype = "log2"
    quanttype = "rlog"
    xval = paste0(quanttype, contrast_base_level, "mean")
    yval = paste0(quanttype, contrast_level_2, "mean")
    l1_lims = c(-1, 11)
    alu_lims = c(-1, 11)
    herv_lims = c(-1,11)

    matches = grep("L1HS:L1", genes)
    l1 = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(l1_lims) + ylim(l1_lims) + 
        scale_color_manual(values = c("red", "blue", "grey")) +
        theme_cowplot() +
        ggtitle("L1HS") +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))

    matches = grep("AluY:Alu", genes)
    alu = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(alu_lims) + ylim(alu_lims) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("AluY") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))


    matches = grep("HERVK", genes)
    herv = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(herv_lims) + ylim(herv_lims) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("HERVK") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))

    agg = results %>% filter(ActiveTE != "Other") %>% ggplot() +
        geom_violin(aes(x = ActiveTE, y = log2FoldChange), draw_quantiles = c(0.5)) +
        geom_hline(aes(yintercept = 0), color = 'red') +
        ggtitle("RTEs") +
        coord_fixed() +
        + ylim(c(-11,11)) +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        theme(aspect.ratio=1)

    legend <- get_legend(
            # create some space to the left of the legend
            l1 + theme(legend.box.margin = margin(0, 0, 0, 1))
            )

    p = plot_grid(agg + theme(legend.position="none"), l1 + theme(legend.position="none"), alu +
            theme(legend.position="none"), herv +
            theme(legend.position="none"), nrow = 1,
            rel_widths = c(1,1,1,1), labels = "AUTO",
            align = "vh",
            axis = 'bt')

    pdf(paste(outputdir, contrast, paste0("activeelementContrastplot", ".pdf"), sep = '/'), width=14, height=4)
    print(plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep=' ')))
    dev.off()

    # options(repr.plot.width = 30, repr.plot.height = 10)
    number_to_sample = 1000
    matches = grep("L1:LINE", genes) %>% sample(number_to_sample,replace=FALSE)
    alpha = 0.9
    l1 = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(l1_lims) + ylim(l1_lims) + 
        scale_color_manual(values = c("red", "blue", "grey")) +
        theme_cowplot() +
        ggtitle("L1") +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))


    matches = grep("Alu:SINE", genes) %>% sample(number_to_sample,replace=FALSE)
    alu = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(alu_lims) + ylim(alu_lims) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("Alu") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))


    matches = grep("ERV.:LTR", genes) %>% sample(number_to_sample,replace=FALSE)
    herv = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(herv_lims) + ylim(herv_lims) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("ERV") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))

    agg = results %>% filter(Family != "Other") %>% ggplot() +
        geom_violin(aes(x = Family, y = log2FoldChange), draw_quantiles = c(0.5)) +
        geom_hline(aes(yintercept = 0), color = 'red') +
        ggtitle("RTEs") +
        coord_fixed() +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) + theme(aspect.ratio=1)



    legend <- get_legend(
            # create some space to the left of the legend
            l1 + theme(legend.box.margin = margin(0, 0, 0, 1))
            )

    p = plot_grid( agg + theme(legend.position="none"),l1 + theme(legend.position="none"), alu +
            theme(legend.position="none"), herv +
            theme(legend.position="none"), nrow = 1,
            rel_widths = c(1,1,1,1), labels = "AUTO",
            align = "vh",
            axis = 'bt')

    pdf(paste(outputdir, contrast, paste0("familyContrastplot", ".pdf"), sep = '/'), width=14, height=4)
    print(plot_grid(p, legend, rel_widths = c(3, .4))+ ggtitle(paste("DE repeats", contrast, sep=' ')))
    dev.off()

    # options(repr.plot.width = 30, repr.plot.height = 10)
    matches = grep("L1HS:L1", genes)

    L1HS = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(l1_lims) + ylim(l1_lims) + 
        scale_color_manual(values = c("red", "blue", "grey")) +
        theme_cowplot() +
        ggtitle("L1HS") +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))


    matches = grep("AluY:Alu", genes)
    aluY = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(alu_lims) + ylim(alu_lims) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("AluY") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))


    matches = grep("HERVK(.)*int", genes)
    hervk = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(herv_lims) + ylim(herv_lims) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("HERVK") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))

    aggactive = results %>% filter(ActiveTE != "Other") %>% ggplot() +
        geom_violin(aes(x = ActiveTE, y = log2FoldChange), draw_quantiles = c(0.5)) +
        geom_hline(aes(yintercept = 0), color = 'red') +
        ggtitle("RTEs") +
        coord_fixed() +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        theme(aspect.ratio=1)


    # options(repr.plot.width = 30, repr.plot.height = 10)


    #######################

    number_to_sample = 1000
    matches = grep("L1:LINE", genes) %>% sample(number_to_sample,replace=FALSE)
    alpha = 0.9
    l1 = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(l1_lims) + ylim(l1_lims) + 
        scale_color_manual(values = c("red", "blue", "grey")) +
        theme_cowplot() +
        ggtitle("L1") +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))


    matches = grep("Alu:SINE", genes) %>% sample(number_to_sample,replace=FALSE)
    alu = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(alu_lims) + ylim(alu_lims) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("Alu") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))


    matches = grep("ERV.:LTR", genes) %>% sample(number_to_sample,replace=FALSE)
    herv = results[matches,] %>% ggplot() +
        geom_point(aes(x = .data[[xval]], y = .data[[yval]], color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(herv_lims) + ylim(herv_lims) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("ERV") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x =paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " ")  , y= paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = ' '))

    agg = results %>% filter(Family != "Other") %>% ggplot() +
        geom_violin(aes(x = Family, y = log2FoldChange), draw_quantiles = c(0.5)) +
        geom_hline(aes(yintercept = 0), color = 'red') +
        ggtitle("RTEs") +
        coord_fixed() +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) + theme(aspect.ratio=1)



    legend <- get_legend(
            # create some space to the left of the legend
            l1 + theme(legend.box.margin = margin(0, 0, 0, 1))
            )

    p = plot_grid( agg + theme(legend.position="none"),l1 + theme(legend.position="none"), alu +
            theme(legend.position="none"), herv +
            theme(legend.position="none"),
            aggactive + theme(legend.position="none"),L1HS + theme(legend.position="none"), aluY +
            theme(legend.position="none"), hervk +
            theme(legend.position="none"),
            nrow = 2,
            rel_widths = c(1,1,1,1), labels = "AUTO",
            align = "vh",
            axis = 'bt')

    pdf(paste(outputdir, contrast, paste0("combinedContrastplot", ".pdf"), sep = '/'), width=14, height=8)
    print(plot_grid(p, legend, rel_widths = c(3, .4))+ ggtitle(paste("DE repeats", contrast, sep=' ')))
    dev.off()

    # options(repr.plot.width = 30, repr.plot.height = 10)

}
}

# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:

# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. Customization of this outlier filtering and description of functionality for replacement of outlier counts and refitting is described below
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. Description and customization of independent filtering is described below


#################### Venn diagrams

for (counttype in snakemake@params[["counttypes"]]) {
    outputdir = paste(snakemake@params[["outputdir"]], counttype, sep = "/")


    UP = list()
    DOWN = list()

    for (contrast in snakemake@params[["contrasts"]]) {
    df = read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast , "results.csv", sep = "/"))
    colnames(df)[1] <- "Geneid"

    dfUP = df[df$pvalue < 0.05 & df$log2FoldChange > 0, 'Geneid'] %>% na.omit()
    dfDOWN = df[df$pvalue < 0.05 & df$log2FoldChange < 0, 'Geneid'] %>% na.omit()

    UP[[gsub("condition_", "", contrast)]] = dfUP
    DOWN[[gsub("condition_", "", contrast)]] = dfDOWN
    }

    fit <- euler(UP, shape = "ellipse")
    upplot = plot(fit, fill = c("tan1", "steelblue1"),
        quantities = TRUE,
        labels = list(font = 4))

    fit <- euler(DOWN, shape = "ellipse")
    downplot = plot(fit, fill = c("tan1", "steelblue1"),
        quantities = TRUE,
        labels = list(font = 4))
        
    plots = list(upplot, downplot)

    grid = plot_grid(plotlist = plots, labels = c("UP", "DOWN"), ncol = 1)
    grid
    pdf(paste(outputdir, paste0("VennDiagram", "allDEGS", ".pdf"), sep = '/'), width=4, height=6)
    print(grid)
    dev.off()

}

###good


save.image()
####################

for (counttype in snakemake@params[["telocaltypes"]]) {
    outputdir = paste(snakemake@params[["outputdir"]], counttype, sep = "/")
    pattern = list("L1:LINE", "Alu:SINE", "ERVK:LTR", "L1HS:L1", "AluY:Alu", "HERVK(.)*int")
    names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int")
    plots = list()
    for (e in names(pattern)) {
        rte = e
        contrastL = list()
        for (contrast in snakemake@params[["contrasts"]]) {
            df = read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast , "results.csv", sep = "/"))
            colnames(df)[1] <- "Geneid"

            matches = grep(pattern[e], df$Geneid)
            searched_element = df[matches,]
            sigsearched = list(searched_element[searched_element$padj < 0.05, 'Geneid'] %>% na.omit())  
            names(sigsearched) <- gsub("condition_", "", contrast)
            contrastL = c(contrastL, sigsearched)
            }
        fit <- euler(contrastL, shape = "ellipse")
        p = plot(fit, fill = c("tan1", "steelblue1"),
            quantities = TRUE,
            labels = list(font = 4),
            main = rte)
        pdf(paste(outputdir, paste0("VennDiagram", rte, ".pdf"), sep = '/'), width=4, height=4)
        print(p)
        dev.off()

        p_no_title = plot(fit, fill = c("tan1", "steelblue1"),
        quantities = TRUE,
        labels = list(font = 4))

        plots = c(plots, list(p_no_title))
    }
grid = plot_grid(plotlist = plots, labels = names(pattern), ncol = 3)
pdf(paste(outputdir, paste0("VennDiagram", "RTEs", ".pdf"), sep = '/'), width=9, height=6)
print(grid)
dev.off()
}




####################
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
geom_histogram(aes(x=length), fill = colors[i]) + theme_cowplot()
plots = c(plots, list(p))
}
head(repmasker)

grid = plot_grid(plotlist = plots, labels = names(pattern))
pdf(paste(outputdir, paste0("histogram", "allRTEs", ".pdf"), sep = '/'), width=10, height=6)
print(grid)
dev.off()


##############

x <- data.frame()
write.table(x, file=snakemake@output[['outfile']], col.names=FALSE)

#############


 