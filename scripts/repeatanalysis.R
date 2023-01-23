log <- file(snakemake@log[[1]], open="wt")
sink(log)

library('readr')
library('dplyr')
library("ggplot2")
library('tibble')
library("RColorBrewer")
library("magrittr")
library("cowplot")
for (contrast in snakemake@params[["contrasts"]]) {
    typesen = ifelse("condition_lsen_vs_prol" == contrast, "late", "early")

    ddsres = read.csv(paste(snakemake@params[["inputdir"]], contrast, "results.csv", sep = "/"))
    ddscounts = read.csv(paste(snakemake@params[["inputdir"]], contrast, "counttable.csv", sep = "/"))
    ddsrlogcounts = read.csv(paste(snakemake@params[["inputdir"]], contrast, "rlogcounts.csv", sep = "/"))
    outputdir = snakemake@params[["outputdir"]]

    colnames(ddsres)[1] <- "Geneid"
    colnames(ddscounts)[1] <- "Geneid"
    colnames(ddsrlogcounts)[1] <- "Geneid"

    ddsrlogcounts = ddsrlogcounts %>% mutate(rlogprolmean = ((SRR6515349 +SRR6515350 +SRR6515351)/3), rlogesenmean = ((SRR6515352 +SRR6515353 +SRR6515354)/3), rloglsenmean = ((SRR6515355 +SRR6515356 +SRR6515357)/3) )
    ddscounts = ddscounts %>% mutate(prolmean = ((SRR6515349 +SRR6515350 +SRR6515351)/3), esenmean = ((SRR6515352 +SRR6515353 +SRR6515354)/3), lsenmean = ((SRR6515355 +SRR6515356 +SRR6515357)/3) )
    avoidzero = 1
    ddscounts = ddscounts %>% mutate(log2prolmean = log2(prolmean + avoidzero), log2esenmean = log2(esenmean + avoidzero), log2lsenmean = log2(lsenmean + avoidzero))
    rlogsub = ddsrlogcounts[,c('Geneid','rlogprolmean','rlogesenmean','rloglsenmean')]
    log2sub = ddscounts[,c('Geneid','prolmean','esenmean','lsenmean','log2prolmean','log2esenmean','log2lsenmean')]
    dflist = list(ddsres, rlogsub,log2sub)
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
    matches = grep("L1HS:L1", genes)

    l1 = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,7)) + ylim(c(-1,7)) + 
        scale_color_manual(values = c("red", "blue", "grey")) +
        theme_cowplot() +
        ggtitle("L1HS") +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))


    matches = grep("AluY:Alu", genes)
    alu = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,7)) + ylim(c(-1,7)) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("AluY") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))


    matches = grep("HERVK", genes)
    herv = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,11)) + ylim(c(-1,11)) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("HERVK") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))

    agg = results %>% ggplot() +
        geom_violin(aes(x = ActiveTE, y = log2FoldChange), draw_quantiles = c(0.5)) +
        geom_hline(aes(yintercept = 0), color = 'red') +
        ggtitle("RTEs") +
        coord_fixed() +
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

    matches = grep("L1:LINE", genes) %>% sample(1500,replace=FALSE)
    alpha = 0.9
    l1 = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,7)) + ylim(c(-1,7)) + 
        scale_color_manual(values = c("red", "blue", "grey")) +
        theme_cowplot() +
        ggtitle("L1") +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))


    matches = grep("Alu:SINE", genes) %>% sample(1000,replace=FALSE)
    alu = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,7)) + ylim(c(-1,7)) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("Alu") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))


    matches = grep("ERV.:LTR", genes) %>% sample(1000,replace=FALSE)
    herv = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,11)) + ylim(c(-1,11)) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("ERV") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))

    agg = results %>% ggplot() +
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
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,7)) + ylim(c(-1,7)) + 
        scale_color_manual(values = c("red", "blue", "grey")) +
        theme_cowplot() +
        ggtitle("L1HS") +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))


    matches = grep("AluY:Alu", genes)
    aluY = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,7)) + ylim(c(-1,7)) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("AluY") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))


    matches = grep("HERVK", genes)
    hervk = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance)) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,11)) + ylim(c(-1,11)) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("HERVK") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))

    aggactive = results %>% ggplot() +
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


    matches = grep("L1:LINE", genes) %>% sample(1500,replace=FALSE)
    alpha = 0.9
    l1 = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,7)) + ylim(c(-1,7)) + 
        scale_color_manual(values = c("red", "blue", "grey")) +
        theme_cowplot() +
        ggtitle("L1") +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))


    matches = grep("Alu:SINE", genes) %>% sample(1000,replace=FALSE)
    alu = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,7)) + ylim(c(-1,7)) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("Alu") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))


    matches = grep("ERV.:LTR", genes) %>% sample(1000,replace=FALSE)
    herv = results[matches,] %>% ggplot() +
        geom_point(aes(x = rlogprolmean, y = rloglsenmean, color = Significance), alpha = alpha) +
        geom_abline(intercept = 0)+ 
        coord_fixed() +
        xlim(c(-1,11)) + ylim(c(-1,11)) +
        scale_color_manual(values = c("red", "blue", "grey")) +
        ggtitle("ERV") +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line=element_blank()) +
        labs(x ="rlog counts proliferating", y= paste("rlog counts", typesen, "senescent", sep = ' '))

    agg = results %>% ggplot() +
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
x <- data.frame()
write.table(x, file=snakemake@output[['outfile']], col.names=FALSE)


# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:

# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. Customization of this outlier filtering and description of functionality for replacement of outlier counts and refitting is described below
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. Description and customization of independent filtering is described below
