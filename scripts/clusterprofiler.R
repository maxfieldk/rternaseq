log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(magrittr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(stringr)
library(cowplot)
library(pathview)
library(ggplot2)
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(ggbeeswarm)

save.image()


basedir = getwd()
for (contrast in snakemake@params[["contrasts"]]) {

    res = read.csv(paste(snakemake@params[["inputdir"]], contrast, "results.csv", sep = "/"))
    outputdir = snakemake@params[["outputdir"]]
    #res = read.csv("/home/mk/scienceL/marco/results/agg/deseq2/condition_lsen_vs_prol/results.csv")
    rownames(res) = res$X
    drop <- c("X")
    res = res[,!(names(res) %in% drop)]

    ## Custom Analysis
    #### Code that applies to both enrichment and GSE:
    senmayo = read.delim(snakemake@params[["SenMayoHuman"]])
    #senmayo = read.delim("/home/mk/scienceL/ref/genesets/SenMayoGeneSetHuman.txt")
    senmayo_mod = mutate(senmayo, setname = "SenMayo")

    #Get all keys
    k <- keys(org.Hs.eg.db,keytype="SYMBOL")
    GOannot = AnnotationDbi::select(org.Hs.eg.db, keys=k, columns=c("GO"), keytype="SYMBOL")
    GOannotshrunk = GOannot[,c("SYMBOL", "GO")]

    customSYMBOL = c(GOannotshrunk$SYMBOL, senmayo_mod$Gene.human.)
    customGO = c(GOannotshrunk$GO, senmayo_mod$setname)
    custom_go2gene = data.frame(go=customGO, gene =customSYMBOL)
    ### Enrichment analysis
    degUP <- rownames(res[((res$padj < 0.05)&(res$log2FoldChange > 0.5)),])
    degDOWN <- rownames(res[((res$padj < 0.05)&(res$log2FoldChange < 0.5)),])
    deg <- rownames(res[((res$padj < 0.05)&(abs(res$log2FoldChange) > 0.5)),])




    top100degs = head(res[order(res$padj),], 100)
    allgenes <- rownames(res)
    oraUP = enricher(degUP, TERM2GENE=custom_go2gene, universe=allgenes)
    oraDOWN = enricher(degDOWN, TERM2GENE=custom_go2gene, universe=allgenes)
    #up
    eg = bitr(oraUP@result$ID, fromType="GOID", toType="TERM", OrgDb="GO.db", drop = FALSE)
    eg[eg$GOID =="SenMayo",]$TERM = "SenMayo"
    term = eg$TERM
    oraUP@result$Description <- term

    #down
    eg = bitr(oraDOWN@result$ID, fromType="GOID", toType="TERM", OrgDb="GO.db", drop = FALSE)
    eg[eg$GOID =="SenMayo",]$TERM = "SenMayo"
    term = eg$TERM
    oraDOWN@result$Description <- term
    oracount = 1
    for (ora in c(oraUP, oraDOWN)) {
        if (oracount == 1) {
        direction = "UP"
        } else {
        direction = "DOWN"
        }
        plot_data = as.data.frame(ora) %>% mutate(qscore = -log(p.adjust, base = 10))
        # Select the columns to use in the plot
        plot_data <- plot_data[, c("Description", "qscore")]

        plot_data$highlight <- ifelse(plot_data$Description == "SenMayo", "highlight", "normal")

        plot_data$Description <- str_wrap(plot_data$Description, 30)

        # Reorder the data by p.adjust
        plot_data <- plot_data[order(-plot_data$qscore), ]
        print(plot_data)

        # Select the top 20 terms
        top_terms1 <- plot_data[1:15, ]
        top_terms2 <- plot_data[16:30, ]

        # Create the barplot
        p1 = ggplot(top_terms1, aes(x = qscore, y = reorder(Description, qscore), fill = highlight)) +
        geom_bar(stat = "identity") +
        labs(x = "-log10 p.adjust", y = "") +
        scale_fill_manual(values = c("normal" = "steelblue", "highlight" = "red")) +
        xlim(c(0, 30)) +
        theme_cowplot(14)

        p2 = ggplot(top_terms2, aes(x = qscore, y = reorder(Description, qscore), fill = highlight)) +
        geom_bar(stat = "identity") +
        labs(x = "-log10 p.adjust", y = "") +
        scale_fill_manual(values = c("normal" = "steelblue", "highlight" = "red")) +
        xlim(c(0, 30)) +
        theme_cowplot(14)


        legend <- get_legend(
        # create some space to the left of the legend
        p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
        )

        prow = plot_grid(p1 + theme(legend.position="none"), p2 +
        theme(legend.position="none"), nrow = 1,
        rel_widths = c(1,1),
        align = 'vh')


        #options(repr.plot.width = 10, repr.plot.height = 6)

        pdf(paste(outputdir, contrast, "hypgeo", paste0("go_enriched", direction, ".pdf"), sep = '/'), width=10, height=6)
        print(plot_grid(prow)+ ggtitle(paste("go enriched",direction, contrast, sep=' ')))
        dev.off()

        oracount = oracount + 1
    }

    ### GSE Analysis
    res <- res[order(-res$stat),]
    ordered_gene_list <- res$stat
    names(ordered_gene_list) <- rownames(res)
    gse = GSEA(ordered_gene_list, TERM2GENE=custom_go2gene, eps = 0)

    eg = bitr(gse@result$ID, fromType="GOID", toType="TERM", OrgDb="GO.db", drop = FALSE)
    eg[eg$GOID =="SenMayo",]$TERM = "SenMayo"
    term = eg$TERM


    gse@result$Description <- term
    gset = "SenMayo"

    pdf(paste(outputdir, contrast, "gsea", paste0("SenMayo", ".pdf"), sep = '/'), width=10, height=6)
    print(gseaplot(gse, geneSetID = gset, title = paste(gse@result[gset,]$Description, contrast, sep = ' ')))
    dev.off()

    pdf(paste(outputdir, contrast, "gsea", paste0("dotplot", ".pdf"), sep = '/'), width=10, height=9)
    print(dotplot(gse, showCategory=30) + ggtitle(paste("GSEA", contrast, sep=' ')) + theme_cowplot())
    dev.off()

    pdf(paste(outputdir, contrast, "gsea", paste0("ridgeplot", ".pdf"), sep = '/'), width=8, height=12)
    print(ridgeplot(gse)+ ggtitle(paste("GSEA", contrast, sep=' ')))
    dev.off()

    gse <- pairwise_termsim(gse)
    pdf(paste(outputdir, contrast, "gsea", paste0("emapplot", ".pdf"), sep = '/'), width=7, height=6)
    print(emapplot(gse, color = "enrichmentScore", showCategory = 20, layout = "nicely")+ ggtitle(paste("GSEA", contrast, sep=' ')))
    dev.off()


    #### KEGG
    eg = bitr(deg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    kk <- enrichKEGG(gene  = eg$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

    pdf(paste(outputdir, contrast, "KEGG", paste0("barplot", ".pdf"), sep = '/'), width=7, height=6)
    print(barplot(kk, showCategory = 20)+ ggtitle(paste("KEGG", 'all DEGs', contrast, sep=' ')))
    dev.off()

    pdf(paste(outputdir, contrast, "KEGG", paste0("dotplot", ".pdf"), sep = '/'), width=7, height=6)
    print(dotplot(kk, showCategory = 20)+ ggtitle(paste("KEGG", 'all DEGs', contrast, sep=' ')))
    dev.off()

    eg = bitr(rownames(res), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    named_foldchange <- setNames(res$log2FoldChange, eg$ENTREZID)
    setwd(paste(outputdir, contrast, "KEGG", "topsigpathways", sep = '/'))
    for (e in head(kk$ID, 20)) {
        pathview(gene.data  = named_foldchange,pathway.id = e,species="hsa", gene.idtype ="entrez", limit = 10)
    }
    setwd(basedir)
    setwd(paste(outputdir, contrast, "KEGG", "importantpathways", sep = '/'))
    my_fav_kegg_pathways = c(tlr="04620", cgas="04623", rlr="04622", senescence="04218")
    for (e in my_fav_kegg_pathways) {
        pathview(gene.data  = named_foldchange,pathway.id = e,species="hsa", gene.idtype ="entrez", limit = 10)
    }

    setwd(basedir)

}

########################################## plotting individual genes
sample_table = read_csv(snakemake@params[["sample_table"]])
#sample_table = read_csv("/users/mkelsey/data/senescence/conf/sample_table.csv")
sample_name = sample_table$sample_name
condition = sample_table$condition
named_list = list()
for (i in seq(length(sample_name))) {
    print(i)
    named_list[[ sample_name[[i]] ]] = condition[[i]]
    i = i+1
}

counts = read_csv(snakemake@input[["normcounttable"]][1])
#counts = read_csv("/users/mkelsey/data/senescence/results/agg/deseq2/star/condition_SEN_vs_PRO/counttablesizenormed.csv")
levels = snakemake@params[["levels"]]
# levels = c(
#       "PRO",
#   "QUI",
#   "SEN",
#   "X3TC",
#   "FTC",
#   "CAS",
#   "KREB"
# )
counts = counts %>% pivot_longer(cols= -1) %>% pivot_wider(names_from = "...1",values_from = "value")
counts = counts %>% plyr::rename(c("name" = "sample"))
counts = counts %>% mutate(condition= unname(unlist(named_list[sample])))

genelists = snakemake@params[["genelistsforplot"]]

for (genelistname in genelists) {
    browser()
    genelist = read_csv(genelistname, col_names = FALSE)
    genelist = genelist$X1
    #genelist = c("IL1A", "IL6", "IFNA1", "IGFBP3", "CGAS")
    plots = list()
    for (gene in genelist) {
        if (gene %in% colnames(counts)) {      
            #first create a grouped summary for the bar plot
            gd <- counts %>% 
                    group_by(condition) %>% 
                    summarise(
                    {{gene}} := mean(.data[[gene]]),
                    )
            #make the plot
            p_no_title = counts %>% ggplot() + 
                geom_col(data = gd, aes(x = factor(condition, level=levels), y = .data[[gene]], fill = factor(condition, level=levels))) + 
                geom_point(aes(x = factor(condition, level=levels), y = .data[[gene]]), position = position_dodge2(width = 0.2)) +
                theme_cowplot() + 
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
                xlab("") + theme(legend.position = "none")
            #pdf(paste(outputdir, "genes", paste0(gene, ".pdf"), sep = '/'), width=4, height=4)
            pdf("zbplot.pdf", width=4, height=4)
            print(p_no_title)
            dev.off()
            #add to plot list
            plots = c(plots, list(p_no_title))
        }

    }

    grid = plot_grid(plotlist = plots, labels = "AUTO", ncol = 3)
    pdf(paste(outputdir, "genes", paste0(genelistname, ".pdf"), sep = '/'), width=9, height=6)
    print(grid)
    dev.off()


    #now make a limited version of plot with fewer categories
    for (contrast in snakemake@params[["contrasts"]]) {
        split = str_split(contrast, "_")
        treatment = split[[1]][2]
        base = split[[1]][4]

        plots_fewer_cat = list()
        for (gene in genelist) {
            if (gene %in% colnames(counts)) {      
                #first create a grouped summary for the bar plot
                gd <- counts %>% 
                    filter(condition %in% c(base, treatment)) %>% group_by(condition) %>%
                    summarise(
                    {{gene}} := mean(.data[[gene]]),
                    )
                p_no_title = counts %>% filter(condition %in% c(base, treatment)) %>% ggplot() + 
                    geom_col(data = gd, aes(x = factor(condition, level=levels), y = .data[[gene]], fill = factor(condition, level=levels))) + 
                    geom_point(aes(x = factor(condition, level=levels), y = .data[[gene]]), position = position_dodge2(width = 0.2)) + 
                    theme_cowplot() + 
                    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
                    xlab("") + 
                    guides(fill = FALSE)
                pdf(paste(outputdir, contrast, "genes", paste0(gene, ".pdf"), sep = '/'), width=3, height=3)
                print(p_no_title)
                dev.off()
                plots_fewer_cat = c(plots_fewer_cat, list(p_no_title))
            }
        }
        grid = plot_grid(plotlist = plots_fewer_cat, labels = "AUTO", ncol = 3)
        pdf(paste(outputdir, contrast, "genes", paste0(genelistname, ".pdf"), sep = '/'), width=9, height=6)
        print(grid)
        dev.off()
    }

}

x <- data.frame()
write.table(x, file=snakemake@output[['outfile']], col.names=FALSE)


