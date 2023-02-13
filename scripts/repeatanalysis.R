# note! there are dataset specific parameters set! need to edit accordingly

log <- file(snakemake@log[[1]], open = "wt")
sink(log)

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
library("ggrepel")
save.image()

peptable <- read.csv(snakemake@params[["peptable"]])
# order matters for the colors!
contrast_colors <- snakemake@params[["contrast_colors"]]
contrast_colors <- unname(unlist(contrast_colors))
condition_colors <- snakemake@params[["condition_colors"]]
condition_colors <- unname(unlist(condition_colors))
repeatsannotated = read_table(snakemake@params[["repeatsannotatted"]], col_names = FALSE) 
telocalrepeatsannotated = read_table(snakemake@params[["telocalmapping"]], col_names = FALSE)
######## MAIN FUNCTIONS
#### GENERAL UTILITY
my_comma <- function(num) {
    format(round(as.numeric(num), 1), nsmall = 0, big.mark = ",")
}
#### PLOTTING
plotRTE <- function(rte, classificationlevel, df, lims = c(-1, 11), title = "", colors = c("red", "orange", "grey"), alpha = 0.9, number_to_sample = NULL, annotate_top = 0) {
    dat <- df %>% filter(.data[[classificationlevel]] == rte)
    if (is.null(number_to_sample)) {
        arranged <- dat %>% arrange(padj)
        arranged$rank <- row.names(arranged)
        arranged$rank <- as.integer(arranged$rank)
        arranged <- arranged %>% mutate(labelshow = ifelse(rank <= annotate_top, "yes", "no"))
        rteplot <- arranged %>% ggplot(aes(x = .data[[xval]], y = .data[[yval]])) +
            geom_point(aes(color = Significance, shape = region2)) +
            geom_abline(intercept = 0) +
            coord_fixed() +
            xlim(lims) +
            ylim(lims) +
            scale_color_manual(values = colors) +
            scale_shape_manual(name = "Genomic\nContext", values = c("Genic" = 19, "Non-Genic" = 17)) +            theme_cowplot() +
            ggtitle(title) +
            panel_border(color = "black", linetype = 1, remove = FALSE) +
            theme(axis.line = element_blank()) +
            labs(x = paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " "), y = paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = " ")) +
            geom_label_repel(
                data = . %>% mutate(label = ifelse(labelshow == "yes", .data[["teorgenename"]], "")),
                aes(label = label),
                box.padding = 0.5,
                point.padding = 0.5,
                nudge_x = -.15,
                nudge_y = 1.5,
                min.segment.length = 0,
                max.overlaps = 2,
                segment.color = "grey50"
            )
        return(rteplot)
    } else {
        number_to_sample <- min(length(rownames(dat)), number_to_sample)
        datasampled <- dat[sample(rownames(dat), number_to_sample, replace = FALSE), ]
        arranged <- datasampled %>% arrange(padj)
        arranged$rank <- row.names(arranged)
        arranged$rank <- as.integer(arranged$rank)
        arranged <- arranged %>% mutate(labelshow = ifelse(rank <= annotate_top, "yes", "no"))
        rteplot <- arranged %>% ggplot(aes(x = .data[[xval]], y = .data[[yval]])) +
            geom_point(aes(color = Significance, shape = region2)) +
            geom_abline(intercept = 0) +
            coord_fixed() +
            xlim(lims) +
            ylim(lims) +
            scale_color_manual(values = colors) +
            scale_shape_manual(name = "Genomic\nContext", values = c("Genic" = 19, "Non-Genic" = 17)) +
            theme_cowplot() +
            ggtitle(title) +
            panel_border(color = "black", linetype = 1, remove = FALSE) +
            theme(axis.line = element_blank()) +
            labs(x = paste(quanttype, "counts", levelslegendmap[contrast_base_level], sep = " "), y = paste(quanttype, "counts", levelslegendmap[contrast_level_2], sep = " ")) +
            geom_label_repel(
                data = . %>% mutate(label = ifelse(labelshow == "yes", .data[["teorgenename"]], "")),
                aes(label = label),
                box.padding = 0.5,
                point.padding = 0.5,
                nudge_x = -.15,
                nudge_y = 1.5,
                min.segment.length = 0,
                max.overlaps = Inf,
                segment.color = "grey50"
            )
        return(rteplot)
    }
}

plotAggRTE <- function(df, classificationlevel, repeattypesallowed = NULL, groupByRegion = TRUE, valuescol = "log2FoldChange") {
    if (is.null(repeattypesallowed)) {
        df <- df %>%
            filter(.data[[classificationlevel]] != "Other")
    } else {
        df <- df %>%
            filter(.data[["Family"]] %in% repeattypesallowed)
    }
    if (groupByRegion) {
    agg <- df %>%
        ggplot(aes(fill = region2, x = .data[[classificationlevel]], y = .data[[valuescol]])) +
        geom_violin(draw_quantiles = c(0.5)) +
        geom_hline(aes(yintercept = 0), color = "red") +
        scale_fill_manual(name = "", values = c("Genic" = "blue", "Non-Genic" = "green")) +
        ggtitle("RTEs") +
        coord_fixed() +
        ylim(c(-11, 11)) +
        theme_cowplot() +
        theme(legend.position = c(0.05, 0.18)) +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line = element_blank()) +
        theme(aspect.ratio = 1)

    } else {
    agg <- df %>%
        ggplot(aes(x = .data[[classificationlevel]], y = .data[[valuescol]])) +
        geom_violin(draw_quantiles = c(0.5)) +
        stat_summary(fun = "mean", geom = "point", color = "black") +
        geom_hline(aes(yintercept = 0), color = "red") +
        ggtitle("RTEs") +
        coord_fixed() +
        ylim(c(-11, 11)) +
        theme_cowplot() +
        panel_border(color = "black", linetype = 1, remove = FALSE) +
        theme(axis.line = element_blank()) +
        theme(aspect.ratio = 1)
    }
    numberxaxis <- df[, classificationlevel] %>%
        unique() %>%
        length()
    if (numberxaxis > 3) {
        agg <- agg + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    }

    return(agg)
}

makeLegendGrob <- function(labels, fill = c("blue", "red"), cex = 1) {
    return(legendGrob(labels,
        pch = 21,
        gp = gpar(
            col = "black",
            fill = fill,
            cex = cex
        )
    ))
}

euPlot <- function(fit, contrasts, fill = factor(c("blue", "red"), ordered = TRUE), alpha = 0.7, legendcex = 1, legend = TRUE, labels = FALSE, main = "") {
    if (legend == TRUE) {
        legd <- makeLegendGrob(contrasts, fill, legendcex)
        p <- plot(
            fit,
            main = main,
            fills = list(fill = fill, alpha = alpha),
            legend = FALSE,
            quantities = list(col = "black"),
            labels = labels
        )
        row <- plot_grid(p)
        plot <- plot_grid(row, legd, nrow = 2, rel_heights = c(1, 0.25))
    } else {
        plot(
            fit,
            main = main,
            fills = list(fill = fill, alpha = alpha),
            legend = FALSE,
            quantities = list(col = "black"),
            labels = labels
        )
    }
}


######
##########################
contrasts <- snakemake@params[["contrasts"]]
levelslegendmap <- snakemake@params[["levelslegendmap"]]
telocaltypes <- snakemake@params[["telocaltypes"]]

for (counttype in telocaltypes) {
    for (contrast in contrasts) {
        ####### DATA WRANGLING
        contrast_level_2 <- unlist(strsplit(contrast, "_", fixed = TRUE))[[2]]
        contrast_base_level <- unlist(strsplit(contrast, "_", fixed = TRUE))[[4]]
        conditions <- c(contrast_base_level, contrast_level_2)

        ddsres <- read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast, "results.csv", sep = "/"))
        ddscounts <- read_csv(paste(snakemake@params[["inputdir"]], counttype, contrast, "counttablesizenormed.csv", sep = "/"))
        ddsrlogcounts <- read_csv(paste(snakemake@params[["inputdir"]], counttype, contrast, "rlogcounts.csv", sep = "/"))

        outputdir <- paste(snakemake@params[["outputdir"]], counttype, sep = "/")

        colnames(ddsres)[1] <- "Geneid"
        colnames(ddscounts)[1] <- "Geneid"
        colnames(ddsrlogcounts)[1] <- "Geneid"




        avoidzero <- 1
        meancols <- c()
        log2meancols <- c()
        rlogmeancols <- c()

        for (condition in conditions) {
            condition_samples <- filter(peptable, condition == {{ condition }})$sample_name
            s <- ""
            for (e in condition_samples) {
                if (s == "") {
                    s <- paste0("`", e, "`")
                } else {
                    s <- paste0(s, "+", "`", e, "`")
                }
            }
            s <- paste0("(", s, ")", "/", length(condition_samples))
            meancol <- paste0(condition, "mean")
            ddscounts <- ddscounts %>% mutate({{ meancol }} := eval(parse(text = s)))
            rlogmeancol <- paste0("rlog", condition, "mean")
            ddsrlogcounts <- ddsrlogcounts %>% mutate({{ rlogmeancol }} := eval(parse(text = s)))
            log2meancol <- paste0("log2", condition, "mean")
            ddscounts <- ddscounts %>% mutate({{ log2meancol }} := log2(.data[[meancol]]))
            meancols <- c(meancols, meancol)
            log2meancols <- c(log2meancols, log2meancol)
            rlogmeancols <- c(rlogmeancols, rlogmeancol)
        }

        columns_to_retain <- c("Geneid", meancols, log2meancols)
        rlogcolumns_to_retain <- c("Geneid", rlogmeancols)

        log2sub <- ddscounts[, columns_to_retain]
        rlogsub <- ddsrlogcounts[, rlogcolumns_to_retain]

        dflist <- list(ddsres, log2sub, rlogsub)

        results <- Reduce(function(x, y) merge(x, y, by = "Geneid"), dflist, accumulate = FALSE)
        results <- results %>% mutate(Significance = ifelse(padj < 0.05, ifelse(padj < 0.001, "< 0.001", "< 0.05"), "> 0.05"))
        results <- results %>% mutate(teorgenename = str_split(Geneid, ":", simplify = TRUE)[, 1])
        colnames(telocalrepeatsannotated) <- c("chr", "start", "stop", "teorgenename", "ignore", "strand", "intronOverlapCount", "exonOverlapCount", "region")
        results = left_join(
        results,
        telocalrepeatsannotated,
        by = "teorgenename",
        copy = FALSE,
        suffix = c(".x", ".y"),
        keep = NULL,
        na_matches = c("na", "never"),
        multiple = "any",
        unmatched = "drop"
        )
        results <- results %>% mutate(region2 = ifelse(region == "exon" | region == "intron", "Genic", "Non-Genic"))
        genes <- results$Geneid


        pattern <- list("L1:LINE", "Alu:SINE", "ERVK:LTR")
        names(pattern) <- c("L1", "Alu", "ERVK")
        classification <- "Family"
        results[classification] <- "Other"
        for (e in names(pattern)) {
            matches <- grep(pattern[[e]], genes)
            results[matches, classification] <- e
        }

        pattern <- list("L1HS:L1", "AluY:Alu", "HERVK(.)*int")
        names(pattern) <- c("L1HS", "AluY", "HERVK-int")
        classification <- "ActiveTE"
        results[classification] <- "Other"
        for (e in names(pattern)) {
            matches <- grep(pattern[[e]], genes)
            results[matches, classification] <- e
        }

        # be sure to have the higher level in ontology first e.g L1:LINE, and lowest level last e.g. L1HS:L1
        # So this will filter out all L1HS from L1:LINE. To get all L1:LINE, you will need to call the family level
        # nevermind, only put same ontology in a column
        pattern <- list("L1HS:L1", "AluY:Alu", "HERVK(.)*int", "L1PA2:L1", "L1PA3:L1", "L1PA4:L1", "L1PA5:L1", "L1PA[6789]{1}:L1")
        names(pattern) <- c("L1HS", "AluY", "HERVK-int", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6-9")
        classification <- "TE"
        results[classification] <- "Other"
        for (e in names(pattern)) {
            matches <- grep(pattern[[e]], genes)
            results[matches, classification] <- e
        }

        ################################################ PLOTING
        ##### plot settings!

        # quanttype = "log2"
        quanttype <- "rlog"
        xval <- paste0(quanttype, contrast_base_level, "mean")
        yval <- paste0(quanttype, contrast_level_2, "mean")
        l1_lims <- c(-1, 11)
        alu_lims <- c(-1, 11)
        herv_lims <- c(-1, 11)

        ##### ACTIVE RTE PLOTS
        # annotated plots
        outputdir <- paste(snakemake@params[["outputdir"]], counttype, sep = "/")

        # scatter plots
        scatterplotsNoAnnotations <- list()
        classificationlevels <- c("TE", "Family")
        for (classificationlevel in classificationlevels) {
            print(classificationlevel)
            types2plot <- results[, classificationlevel] %>% unique()
            for (e in types2plot) {
                print(e)
                if (e != "Other") {
                    tempplot <- plotRTE(e, classificationlevel, results, lims = l1_lims, title = e, number_to_sample = 1000, annotate_top = 0)
                    dirname <- file.path(outputdir, contrast, e)	
                    dir.create(dirname, recursive = TRUE)                           
                    pdf(paste(dirname, paste0(e, "NotAnnotated", ".pdf"), sep = "/"), width = 5, height = 4)
                    print(tempplot)
                    dev.off()
                    scatterplotsNoAnnotations[[e]] <- tempplot
                }
            }
        }

        scatterplotsHighlyAnnotated <- list()
        classificationlevels <- c("TE", "Family")
        for (classificationlevel in classificationlevels) {
            print(classificationlevel)
            types2plot <- results[, classificationlevel] %>% unique()
            for (e in types2plot) {
                print(e)
                if (e != "Other") {                    
                    tempplot <- plotRTE(e, classificationlevel, results, lims = l1_lims, title = e, number_to_sample = 1000, annotate_top = 5)
                    dirname <- file.path(outputdir, contrast, e)	
                    dir.create(dirname, recursive = TRUE)                           
                    pdf(paste(dirname, paste0(e, "highlyannotated", ".pdf"), sep = "/"), width = 5, height = 4)
                    print(tempplot)
                    dev.off()
                    scatterplotsHighlyAnnotated[[e]] <- tempplot
                }
            }
        }

        scatterplots <- list()
        classificationlevels <- c("TE", "Family")
        for (classificationlevel in classificationlevels) {
            print(classificationlevel)
            types2plot <- results[, classificationlevel] %>% unique()
            for (e in types2plot) {
                print(e)
                if (e != "Other") {
                    tempplot <- plotRTE(e, classificationlevel, results, lims = l1_lims, title = e, number_to_sample = 1000, annotate_top = 3)
                    dirname <- file.path(outputdir, contrast, e)	
                    dir.create(dirname, recursive = TRUE)                    
                    pdf(paste(dirname, paste0(e, "annotated", ".pdf"), sep = "/"), width = 5, height = 4)
                    print(tempplot)
                    dev.off()
                    scatterplots[[e]] <- tempplot
                }
            }
        }

        # violin plots
        violinplots <- list()
        classificationlevels <- c("TE", "Family", "ActiveTE")
        for (classificationlevel in classificationlevels) {
            tempplot <- plotAggRTE(results, classificationlevel, groupByRegion = TRUE)
            dirname <- file.path(outputdir, contrast)	
            dir.create(dirname, recursive = TRUE)            
            pdf(paste(dirname, paste0("Agg", classificationlevel, ".pdf"), sep = "/"), width = 5, height = 4)
            print(tempplot)
            dev.off()
            violinplots[[classificationlevel]] <- tempplot
        }

        L1violinplots <- list()
        classificationlevels <- c("TE")
        for (classificationlevel in classificationlevels) {
            tempplot <- plotAggRTE(results, classificationlevel, repeattypesallowed = c("L1"), groupByRegion = TRUE)
            dirname <- file.path(outputdir, contrast)	
            dir.create(dirname, recursive = TRUE)            
            pdf(paste(dirname, paste0("AggL1", classificationlevel, ".pdf"), sep = "/"), width = 5, height = 4)
            print(tempplot)
            dev.off()
            L1violinplots[[classificationlevel]] <- tempplot
        }

        ####### region plots

        regionplots <- list()
        classificationlevels <- c("TE", "Family", "ActiveTE")
        for (classificationlevel in classificationlevels) {
            tempplot <- plotAggRTE(results, classificationlevel)
            dirname <- file.path(outputdir, contrast)	
            dir.create(dirname, recursive = TRUE)            
            pdf(paste(dirname, paste0("Agg", classificationlevel, ".pdf"), sep = "/"), width = 5, height = 4)
            print(tempplot)
            dev.off()
            violinplots[[classificationlevel]] <- tempplot
        }


        ########### cow plots
        ###### First with annotations
        # Active RTE
        legend <- get_legend(scatterplots[["L1HS"]] + theme(legend.box.margin = margin(0, 0, 0, 1)))
        p <- plot_grid(violinplots[["ActiveTE"]],
            scatterplots[["L1HS"]] + theme(legend.position = "none"),
            scatterplots[["AluY"]] + theme(legend.position = "none"),
            scatterplots[["HERVK-int"]] + theme(legend.position = "none"),
            nrow = 1,
            rel_widths = c(1, 1, 1, 1), labels = "AUTO",
            align = "vh",
            axis = "bt"
        )
        p <- plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep = " "))
        dirname <- file.path(outputdir, contrast)	
        dir.create(dirname, recursive = TRUE)        
        pdf(paste(dirname, paste0("activeelementContrastplot", ".pdf"), sep = "/"), width = 14, height = 4)
        print(p)
        dev.off()

        # Family
        legend <- get_legend(scatterplots[["L1"]] + theme(legend.box.margin = margin(0, 0, 0, 1)))
        p <- plot_grid(violinplots[["Family"]],
            scatterplots[["L1"]] + theme(legend.position = "none"),
            scatterplots[["Alu"]] + theme(legend.position = "none"),
            scatterplots[["ERVK"]] + theme(legend.position = "none"),
            nrow = 1,
            rel_widths = c(1, 1, 1, 1), labels = "AUTO",
            align = "vh",
            axis = "bt"
        )
        p <- plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep = " "))
        dirname <- file.path(outputdir, contrast)	
        dir.create(dirname, recursive = TRUE)
        pdf(paste(dirname, paste0("FamilyContrastplot", ".pdf"), sep = "/"), width = 14, height = 4)
        print(p)
        dev.off()


        ########## Plot activete and family together
        legend <- get_legend(scatterplots[["L1"]] + theme(legend.box.margin = margin(0, 0, 0, 1)))
        p <- plot_grid(violinplots[["Family"]],
            scatterplots[["L1"]] + theme(legend.position = "none"),
            scatterplots[["Alu"]] + theme(legend.position = "none"),
            scatterplots[["ERVK"]] + theme(legend.position = "none"),
            violinplots[["ActiveTE"]],
            scatterplots[["L1HS"]] + theme(legend.position = "none"),
            scatterplots[["AluY"]] + theme(legend.position = "none"),
            scatterplots[["HERVK-int"]] + theme(legend.position = "none"),
            nrow = 2,
            rel_widths = c(1, 1, 1, 1), labels = "AUTO",
            align = "vh",
            axis = "bt"
        )
        p <- plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep = " "))
        dirname <- file.path(outputdir, contrast)	
        dir.create(dirname, recursive = TRUE)        
        pdf(paste(dirname, paste0("CombinedContrastPlot", ".pdf"), sep = "/"), width = 14, height = 8)
        print(p)
        dev.off()

        ########## L1 together plot
        legend <- get_legend(scatterplots[["L1"]] + theme(legend.box.margin = margin(0, 0, 0, 1)))
        p <- plot_grid(L1violinplots[["TE"]],
            scatterplots[["L1"]] + theme(legend.position = "none"),
            scatterplots[["L1HS"]] + theme(legend.position = "none"),
            scatterplots[["L1PA2"]] + theme(legend.position = "none"),
            scatterplots[["L1PA3"]] + theme(legend.position = "none"),
            scatterplots[["L1PA4"]] + theme(legend.position = "none"),
            scatterplots[["L1PA5"]] + theme(legend.position = "none"),
            scatterplots[["L1PA6-9"]] + theme(legend.position = "none"),
            nrow = 2,
            rel_widths = c(1, 1, 1, 1), labels = "AUTO",
            align = "vh",
            axis = "bt"
        )
        p <- plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep = " "))
        dirname <- file.path(outputdir, contrast)	
        dir.create(dirname, recursive = TRUE)        
        pdf(paste(dirname, paste0("CombinedContrastPlotL1", ".pdf"), sep = "/"), width = 14, height = 8)
        print(p)
        dev.off()


        ############## Same together plots but without the annotations
        # Active RTE
        legend <- get_legend(scatterplotsNoAnnotations[["L1HS"]] + theme(legend.box.margin = margin(0, 0, 0, 1)))
        p <- plot_grid(violinplots[["ActiveTE"]],
            scatterplotsNoAnnotations[["L1HS"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["AluY"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["HERVK-int"]] + theme(legend.position = "none"),
            nrow = 1,
            rel_widths = c(1, 1, 1, 1), labels = "AUTO",
            align = "vh",
            axis = "bt"
        )
        p <- plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep = " "))
        dirname <- file.path(outputdir, contrast)	
        dir.create(dirname, recursive = TRUE)
        pdf(paste(dirname, paste0("activeelementContrastplotNoAnnotations", ".pdf"), sep = "/"), width = 14, height = 4)
        print(p)
        dev.off()

        # Family
        legend <- get_legend(scatterplotsNoAnnotations[["L1"]] + theme(legend.box.margin = margin(0, 0, 0, 1)))
        p <- plot_grid(violinplots[["Family"]],
            scatterplotsNoAnnotations[["L1"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["Alu"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["ERVK"]] + theme(legend.position = "none"),
            nrow = 1,
            rel_widths = c(1, 1, 1, 1), labels = "AUTO",
            align = "vh",
            axis = "bt"
        )
        p <- plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep = " "))
        dirname <- file.path(outputdir, contrast)	
        dir.create(dirname, recursive = TRUE)
        pdf(paste(dirname, paste0("FamilyContrastplotNoAnnotations", ".pdf"), sep = "/"), width = 14, height = 4)
        print(p)
        dev.off()


        ########## Plot activete and family together
        legend <- get_legend(scatterplotsNoAnnotations[["L1"]] + theme(legend.box.margin = margin(0, 0, 0, 1)))
        p <- plot_grid(violinplots[["Family"]],
            scatterplotsNoAnnotations[["L1"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["Alu"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["ERVK"]] + theme(legend.position = "none"),
            violinplots[["ActiveTE"]],
            scatterplotsNoAnnotations[["L1HS"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["AluY"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["HERVK-int"]] + theme(legend.position = "none"),
            nrow = 2,
            rel_widths = c(1, 1, 1, 1), labels = "AUTO",
            align = "vh",
            axis = "bt"
        )
        p <- plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep = " "))
        dirname <- file.path(outputdir, contrast)	
        dir.create(dirname, recursive = TRUE)
        pdf(paste(dirname, paste0("CombinedContrastPlotNoAnnotations", ".pdf"), sep = "/"), width = 14, height = 8)
        print(p)
        dev.off()

        ########## L1 together plot
        legend <- get_legend(scatterplotsNoAnnotations[["L1"]] + theme(legend.box.margin = margin(0, 0, 0, 1)))
        p <- plot_grid(L1violinplots[["TE"]],
            scatterplotsNoAnnotations[["L1"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["L1HS"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["L1PA2"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["L1PA3"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["L1PA4"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["L1PA5"]] + theme(legend.position = "none"),
            scatterplotsNoAnnotations[["L1PA6-9"]] + theme(legend.position = "none"),
            nrow = 2,
            rel_widths = c(1, 1, 1, 1), labels = "AUTO",
            align = "vh",
            axis = "bt"
        )
        p <- plot_grid(p, legend, rel_widths = c(3, .4)) + ggtitle(paste("DE repeats", contrast, sep = " "))
        pdf(paste(outputdir, contrast, paste0("CombinedContrastPlotL1NoAnnotations", ".pdf"), sep = "/"), width = 14, height = 8)
        print(p)
        dev.off()


        ################################################ END#PLOTING
    }
}

# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:

# If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. Customization of this outlier filtering and description of functionality for replacement of outlier counts and refitting is described below
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. Description and customization of independent filtering is described below


################################################ Venn diagrams of non-RTEs
for (counttype in snakemake@params[["counttypes"]]) {
    outputdir <- paste(snakemake@params[["outputdir"]], counttype, sep = "/")
    UP <- list()
    DOWN <- list()

    eulerr_options(quantities = list(font = 2, cex = 1.5), legend = list(side = "bottom", cex = 1.5))

    for (contrast in contrasts) {
        df <- read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast, "results.csv", sep = "/"))
        colnames(df)[1] <- "Geneid"
        dfUP <- df[df$pvalue < 0.05 & df$log2FoldChange > 0, "Geneid"] %>% na.omit()
        dfDOWN <- df[df$pvalue < 0.05 & df$log2FoldChange < 0, "Geneid"] %>% na.omit()
        UP[[gsub("condition_", "", contrast)]] <- dfUP
        DOWN[[gsub("condition_", "", contrast)]] <- dfDOWN
    }
    fit <- euler(UP, shape = "ellipse")
    upplot <- euPlot(fit, gsub("condition_", "", contrasts), fill = contrast_colors, main = "UP", legend = FALSE)

    fit <- euler(DOWN, shape = "ellipse")
    downplot <- euPlot(fit, gsub("condition_", "", contrasts), fill = contrast_colors, main = "DOWN", legend = FALSE)
    legd <- makeLegendGrob(gsub("condition_", "", contrasts))

    plots <- list(upplot, downplot)
    grid <- plot_grid(plotlist = plots, labels = "AUTO", ncol = 1)
    fgrid <- plot_grid(grid, legd, nrow = 2, rel_heights = c(1, 0.2))
    pdf(paste(outputdir, paste0("VennDiagram", "allDEGS", ".pdf"), sep = "/"), width = 4, height = 6)
    print(fgrid)
    dev.off()
}

################################################
master_deRTEs <- list()
for (direction in c("UP", "DOWN")) {
    counttype_deRTEs <- list()
    for (counttype in snakemake@params[["telocaltypes"]]) {
        outputdir <- paste(snakemake@params[["outputdir"]], counttype, sep = "/")
        pattern <- list("L1:LINE", "Alu:SINE", "ERVK:LTR", "L1HS:L1", "AluY:Alu", "HERVK(.)*int", "L1PA(.){0,2}:L1", "L1PA2:L1", "L1PA3:L1", "L1PA4:L1", "L1PA5:L1", "L1PA[6789]{1}:L1")
        names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int", "L1PA", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6-9")
        deRTEs <- list()
        for (e in names(pattern)) {
            rte <- e
            contrastL <- list()
            for (contrast in contrasts) {
                df <- read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast, "results.csv", sep = "/"))
                colnames(df)[1] <- "Geneid"
                matches <- grep(pattern[e], df$Geneid)
                searched_element <- df[matches, ]
                if (direction == "UP") {
                    sigsearched <- list(searched_element[searched_element$padj < 0.05 & searched_element$log2FoldChange > 0, "Geneid"] %>% na.omit())
                } else {
                    sigsearched <- list(searched_element[searched_element$padj < 0.05 & searched_element$log2FoldChange < 0, "Geneid"] %>% na.omit())
                }
                names(sigsearched) <- gsub("condition_", "", contrast)
                contrastL <- c(contrastL, sigsearched)
            }
            deRTEs[[e]] <- contrastL
        }
        counttype_deRTEs[[counttype]] <- deRTEs
    }
    master_deRTEs[[direction]] <- counttype_deRTEs
}


####### Write the DE results to table
mapper <- read.delim(snakemake@params[["telocalmapping"]], sep = "\t", header = FALSE)
#### functions
getPos <- function(rte, tab) {
    rterow <- mapper %>%
        filter(V4 == rte) %>%
        group_by(V4) %>%
        summarize(
            chr = dplyr::first(V1),
            start = min(V2),
            stop = max(V3),
            name = dplyr::first(V4),
            bs = dplyr::first(V5),
            strand = dplyr::first(V6),
            intron = max(V7),
            exon = max(V8)
        )
    return(list(chr = rterow$chr, start = as.numeric(rterow$start), stop = as.numeric(rterow$stop), strand = rterow$strand, intron = as.numeric(rterow$intron), exon = as.numeric(rterow$exon)))
}


dedf <- data.frame(tetype = character(), te = character(), chr = character(), start = numeric(), stop = numeric(), strand = character(), direction = character(), contrast = character(), counttype = character(), length = numeric(), intron = character(), exon = character(), region = character(), stringsAsFactors = FALSE)
i <- 1
rtenames <- c("L1HS", "AluY", "HERVK-int", "L1PA2", "L1PA3", "L1PA4")
for (direction in c("UP", "DOWN")) {
    for (name in rtenames) {
        for (telocaltype in telocaltypes) {
            for (contrast in contrasts) {
                elements <- master_deRTEs[[direction]][[telocaltype]][[name]][[gsub("condition_", "", contrast)]]
                for (rte in elements) {
                    rte <- str_split(rte, ":")[[1]][1]
                    print(rte)
                    location <- getPos(rte, tab)
                    if (is.null(location)) {
                        next
                    }
                    chr <- location["chr"][[1]]
                    start <- location["start"][[1]]
                    stop <- location["stop"][[1]]
                    strand <- location["strand"][[1]]
                    length <- stop - start
                    intron <- location["intron"][[1]]
                    exon <- location["exon"][[1]]
                    region <- ifelse(intron > 0, "intron", ifelse(exon > 0, "exon", "other"))
                    vec <- c(name, rte, chr, start, stop, strand, direction, contrast, telocaltype, length, intron, exon, region)
                    dedf[i, ] <- vec
                    i <- i + 1
                }
            }
        }
    }
}

write.table(dedf, file = snakemake@output[["DETEsbyContrast"]], quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

head(dedf)
outputdir <- paste(snakemake@params[["outputdir"]], sep = "/")

for (counttype in telocaltypes) {
    for (contrast in contrasts) {
        df <- dedf[dedf$counttype == counttype & dedf$contrast == contrast, ]
        pl <- df %>% ggplot() +
            geom_bar(aes(x = region, fill = region)) +
            ggtitle("DE RTE") +
            theme(aspect.ratio = 1) +
            theme_cowplot()
        dirname <- file.path(outputdir, counttype, contrast)
        basename <- "derteByLocation.pdf"
        pdf(file.path(dirname, basename), width = 4, height = 3)
        print(pl)
        dev.off()
        for (tetype in df$tetype %>% unique()) {
            for (direction in df$direction %>% unique()) {
                df2 <- df[df$tetype == tetype & df$direction == direction, ]
                pl <- df2 %>% ggplot() +
                    geom_bar(aes(x = region, fill = region)) +
                    ggtitle(paste("DE", tetype, direction)) +
                    theme(aspect.ratio = 1) +
                    theme_cowplot()
                dirname <- file.path(outputdir, counttype, contrast, tetype)	
                dir.create(dirname, recursive = TRUE)
                basename <- paste0("derteByLocation", tetype, direction, ".pdf")
                pdf(file.path(dirname, basename), width = 4, height = 3)
                print(pl)
                dev.off()
            }
        }
    }
}



########### Are de RTEs the same?

# this will be: for each counttype, what RTEs were DE in all contrasts
shared_des_master <- list()
for (direction in c("UP", "DOWN")) {
    des_by_counttype <- list()
    for (counttype in snakemake@params[["telocaltypes"]]) {
        outputdir <- paste(snakemake@params[["outputdir"]], counttype, sep = "/")
        pattern <- list("L1:LINE", "Alu:SINE", "ERVK:LTR", "L1HS:L1", "AluY:Alu", "HERVK(.)*int", "L1PA(.){0,2}:L1", "L1PA2:L1", "L1PA3:L1", "L1PA4:L1", "L1PA5:L1", "L1PA[6789]{1}:L1")
        names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int", "L1PA", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6-9")
        des <- list()
        for (e in names(pattern)) {
            decontrast <- list()
            for (contrast in contrasts) {
                decontrast <- c(decontrast, list(master_deRTEs[[direction]][[counttype]][[e]][[gsub("condition_", "", contrast)]]))
            }
            set_des <- Reduce(intersect, decontrast)
            des[[e]] <- set_des
        }
        des_by_counttype[[counttype]] <- des
    }
    shared_des_master[[direction]] <- des_by_counttype
}

# this will be: for each RTE, which were DE in all contrasts AND for all count methods
des_by_rte_master <- list()
for (direction in c("UP", "DOWN")) {
    des_by_rte <- list()
    for (e in names(pattern)) {
        des <- list()
        for (countype in snakemake@params[["telocaltypes"]]) {
            des <- c(des, list(shared_des_master[[direction]][[counttype]][[e]]))
        }
        set_des <- Reduce(intersect, des)
        des_by_rte[[e]] <- set_des
    }
    des_by_rte_master[[direction]] <- des_by_rte
}

############ Calculating significance of overalps

## For stats purposes
telocalcountsample1 <- read.csv(snakemake@input[["telocal"]][[1]])
# telocalcountsample1 =read.delim("/users/mkelsey/data/marco/outs/SRR6515349/TElocal/SRR6515349_uniq.cntTable", sep = "\t")
allgenes <- telocalcountsample1[[1]]
pattern <- list("L1:LINE", "Alu:SINE", "ERVK:LTR", "L1HS:L1", "AluY:Alu", "HERVK(.)*int", "L1PA(.){0,2}:L1", "L1PA2:L1", "L1PA3:L1", "L1PA4:L1", "L1PA5:L1", "L1PA[6789]{1}:L1")
names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int", "L1PA", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6-9")

numberOfInstances <- list()
for (e in names(pattern)) {
    rte <- e
    matches <- grep(pattern[e], allgenes)
    numberOfInstances[e] <- length(matches)
}

sigdf <- data.frame(counttype = character(), element = character(), sig = numeric(), universe = character(), direction = character())
i <- 1

for (direction in c("UP", "DOWN")) {
    for (telocaltype in telocaltypes) {
        for (e in names(pattern)) {
            ntotal <- unname(unlist(numberOfInstances[e]))
            contrast_lengths <- c()
            for (thing in master_deRTEs[[direction]][[telocaltype]][[e]]) {
                contrast_lengths <- c(contrast_lengths, length(thing))
            }
            overlap <- length(shared_des_master[[direction]][[telocaltype]][[e]])
            if (length(contrasts) == 2) {
                significance <- phyper(overlap - 1, contrast_lengths[2], ntotal - contrast_lengths[2], contrast_lengths[1], lower.tail = FALSE)
            } else {
                significance <- 75
            }

            vec <- list(telocaltype, e, significance, ntotal, direction)
            sigdf[i, ] <- vec
            i <- i + 1
        }
    }
}
sigdf <- sigdf %>% mutate(sigcode = ifelse(sig == 75, " ", ifelse(sig > 0.05, "NS", ifelse(sig > 0.01, "*", ifelse(sig > 0.001, "**", "***")))))


###########################
for (direction in c("UP", "DOWN")) {
    for (counttype in snakemake@params[["telocaltypes"]]) {
        outputdir <- paste(snakemake@params[["outputdir"]], counttype, sep = "/")
        pattern <- list("L1:LINE", "Alu:SINE", "ERVK:LTR", "L1HS:L1", "AluY:Alu", "HERVK(.)*int", "L1PA(.){0,2}:L1", "L1PA2:L1", "L1PA3:L1", "L1PA4:L1", "L1PA5:L1", "L1PA[6789]{1}:L1")
        names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int", "L1PA", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6-9")
        plots <- list()
        for (e in names(pattern)) {
            rte <- e
            contrastL <- list()
            for (contrast in contrasts) {
                df <- read.csv(paste(snakemake@params[["inputdir"]], counttype, contrast, "results.csv", sep = "/"))
                colnames(df)[1] <- "Geneid"
                matches <- grep(pattern[e], df$Geneid)
                searched_element <- df[matches, ]
                if (direction == "UP") {
                    sigsearched <- list(searched_element[searched_element$padj < 0.05 & searched_element$log2FoldChange > 0, "Geneid"] %>% na.omit())
                } else {
                    sigsearched <- list(searched_element[searched_element$padj < 0.05 & searched_element$log2FoldChange < 0, "Geneid"] %>% na.omit())
                }
                names(sigsearched) <- gsub("condition_", "", contrast)
                contrastL <- c(contrastL, sigsearched)
            }
            fit <- euler(contrastL, shape = "ellipse")
            p <- euPlot(fit, gsub("condition_", "", contrasts), fill = contrast_colors, main = paste(rte, direction, sep = " "), legendcex = 1)
            sigcode <- sigdf[sigdf$counttype == counttype & sigdf$element == e & sigdf$direction == direction, ]$sigcode
            if (sigcode == "NS") {
                size <- 4
            } else {
                size <- 15
            }
            universe <- sigdf[sigdf$counttype == counttype & sigdf$element == e & sigdf$direction == direction, ]$universe
            pl <- ggdraw(p) +
                annotate("text", label = sigcode, x = Inf, y = Inf, vjust = 1, hjust = 1, size = size) +
                annotate("text", label = paste0("\nU: ", my_comma(universe)), x = Inf, y = Inf, vjust = 1, hjust = 1, size = 4)
            dirname <- file.path(outputdir, counttype, tetype)	
            dir.create(dirname, recursive = TRUE)
            pdf(paste(dirname, paste0("VennDiagram", rte, direction, ".pdf"), sep = "/"), width = 4, height = 4)
            print(pl)
            dev.off()
            p_no_title <- euPlot(fit, contrasts, main = rte, fill = contrast_colors, legend = FALSE)
            plnotitle <- ggdraw(p_no_title) +
                annotate("text", label = sigcode, x = Inf, y = Inf, vjust = 1, hjust = 1, size = size) +
                annotate("text", label = paste0("\nU: ", my_comma(universe)), x = Inf, y = Inf, vjust = 1, hjust = 1, size = 4)
            plots[[e]] <- plnotitle
        }

        grid <- plot_grid(plotlist = plots[1:6], ncol = 3)
        legd <- makeLegendGrob(gsub("condition_", "", contrasts), fill = contrast_colors)
        title <- ggdraw() + draw_label(paste("DE", direction, sep = " "), fontface = "bold", size = 20)
        fgrid <- plot_grid(title, grid, legd, nrow = 3, rel_heights = c(0.15, 1, 0.2))
        dirname <- file.path(outputdir, counttype)	
        dir.create(dirname, recursive = TRUE)
        pdf(paste(dirname, paste0("VennDiagram", "RTEs", direction, ".pdf"), sep = "/"), width = 9, height = 6)
        print(fgrid)
        dev.off()

        grid <- plot_grid(plotlist = plots[7:12], ncol = 3)
        legd <- makeLegendGrob(gsub("condition_", "", contrasts), fill = contrast_colors)
        title <- ggdraw() + draw_label(paste("DE", direction, sep = " "), fontface = "bold", size = 20)
        fgrid <- plot_grid(title, grid, legd, nrow = 3, rel_heights = c(0.15, 1, 0.2))
        dirname <- file.path(outputdir)	
        dir.create(dirname, recursive = TRUE)
        pdf(paste(dirname, paste0("VennDiagram", "L1PAs", direction, ".pdf"), sep = "/"), width = 9, height = 6)
        print(fgrid)
        dev.off()

        grid <- plot_grid(plotlist = plots[1:12], ncol = 3)
        legd <- makeLegendGrob(gsub("condition_", "", contrasts), fill = contrast_colors)
        title <- ggdraw() + draw_label(paste("DE", direction, sep = " "), fontface = "bold", size = 20)
        fgrid <- plot_grid(title, grid, legd, nrow = 3, rel_heights = c(0.15, 2, 0.2))
        dirname <- file.path(outputdir)	
        dir.create(dirname, recursive = TRUE)
        pdf(paste(dirname, paste0("VennDiagram", "AllRTEs", direction, ".pdf"), sep = "/"), width = 9, height = 12)
        print(fgrid)
        dev.off()
    }
}


################# write table for DEs shared amongst all constrasts

dedf <- data.frame(tetype = character(), te = character(), chr = character(), start = numeric(), stop = numeric(), strand = character(), direction = character(), length = numeric(), stringsAsFactors = FALSE)
i <- 1
rtenames <- c("L1HS", "AluY", "HERVK-int")
for (direction in c("UP", "DOWN")) {
    for (name in rtenames) {
        dertelist <- des_by_rte_master[[direction]][[name]]
        print(name)
        for (rte in dertelist) {
            rte <- str_split(rte, ":")[[1]][1]
            print(rte)
            location <- getPos(rte, tab)
            if (is.null(location)) {
                next
            }
            chr <- location["chr"][[1]]
            start <- location["start"][[1]]
            stop <- location["stop"][[1]]
            strand <- location["strand"][[1]]
            length <- stop - start
            vec <- c(name, rte, chr, start, stop, strand, direction, length)
            dedf[i, ] <- vec
            i <- i + 1
        }
    }
}
write.table(dedf, file = snakemake@output[["sharedamongallcontrasts_derte"]], quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

#################### RTE in genome histograms
outputdir <- snakemake@params[["outputdir"]]
repeats <- snakemake@params[["repeats"]]
repmasker <- read_delim(repeats, col_names = FALSE) %>%
    mutate_at(c("X2", "X3"), as.numeric) %>%
    mutate(length = X3 - X2)


pattern <- c("LINE/L1", "SINE/Alu", "LTR/ERVK", "L1HS", "AluY", "HERVK(.)*int")
names(pattern) <- c("L1", "Alu", "ERVK", "L1HS", "AluY", "HERVK-int")
activeLengthReq <- c(1, 0, 0, 6000, 290, 7000)
maxlength <- c(7000, 400, 15000, 7000, 400, 15000)
colors <- c("deeppink4", "deepskyblue4", "darkorange4", "deeppink1", "deepskyblue1", "darkorange1")

plots <- list()

for (i in seq(length(pattern))) {
    print(i)
    print(pattern[i])
    p <- repmasker %>%
        filter(length < maxlength[i]) %>%
        filter(grepl(pattern[i], X4)) %>%
        ggplot() +
        geom_histogram(aes(x = length), fill = colors[i]) +
        ggtitle(names(pattern)[i]) +
        theme_cowplot()
    plots <- c(plots, list(p))
}

grid <- plot_grid(plotlist = plots, labels = "AUTO")
pdf(paste(outputdir, paste0("histogram", "allRTEs", ".pdf"), sep = "/"), width = 10, height = 6)
print(grid)
dev.off()

##############

x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)
save.image()
#############
