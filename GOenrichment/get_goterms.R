library(biomaRt)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(GO.db)

Mgenes <- read.table("276genes.txt", header = FALSE, stringsAsFactors = FALSE)$V1
fungi_mart <- useMart("fungi_mart", host = "https://fungi.ensembl.org", dataset="pgraminis_eg_gene")

listAttributes(fungi_mart)
results <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'go_id', 'namespace_1003', 'name_1006', 'definition_1006'),
                 filters = 'ensembl_gene_id',
                 values = genes,
                 mart = fungi_mart)

results_df <- results |>
  filter(go_id != "") |>
  separate_rows(go_id, sep = ";") |>
  dplyr::select(ensembl_gene_id, go_id)

write.table(results_df, file = "output.txt", sep = "\t", row.names = FALSE)

xls_results <- results |>
  filter(ensembl_gene_id %in% Mgenes) |>
  filter(go_id != "") |>
  filter(namespace_1003 %in% c('biological_process','molecular_function','cellular_component'))

modified_results <- xls_results |>
  mutate(prefix = case_when(
    namespace_1003 == "biological_process" ~ "BP",
    namespace_1003 == "molecular_function" ~ "MF",
    namespace_1003 == "cellular_component" ~ "CC"
  )) |>
  mutate(
    go_id = paste(prefix, go_id, sep = ":"),
    name_1006 = paste(prefix, name_1006, sep = ":")
  )

combined_results <- modified_results |>
  group_by(ensembl_gene_id) |>
  summarize(
    chromosome_name = dplyr::first(chromosome_name),
    start_position = dplyr::first(start_position),
    end_position = dplyr::first(end_position),
    go_id = paste(unique(go_id), collapse = "; "),
    name_1006 = paste(unique(name_1006), collapse = "; ")
  )

write.table(combined_results, file="Mgenes_annotations.tsv",sep="\t",row.names=FALSE)

all_xls_results <- results |> filter(ensembl_gene_id %in% Mgenes)
d1 <- unique(all_xls_results$ensembl_gene_id)
Mgenes[!(Mgenes %in% d1)]

all_modified_results <- all_xls_results |>
  mutate(prefix = case_when(
    namespace_1003 == "biological_process" ~ "BP",
    namespace_1003 == "molecular_function" ~ "MF",
    namespace_1003 == "cellular_component" ~ "CC"
  )) |>
  mutate(
    go_id = paste(prefix, go_id, sep = ":"),
    name_1006 = paste(prefix, name_1006, sep = ":")
  )

all_combined_results <- all_modified_results |>
  group_by(ensembl_gene_id) |>
  summarize(
    chromosome_name = dplyr::first(chromosome_name),
    start_position = dplyr::first(start_position),
    end_position = dplyr::first(end_position),
    go_id = paste(unique(go_id), collapse = "; "),
    name_1006 = paste(unique(name_1006), collapse = "; ")
  )

write.table(all_combined_results, file="all_Mgenes_annotations.tsv",sep="\t",row.names=FALSE)

################################################################################
library(AnnotationHub)
library(shiny)

ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")

orgdb <- query(orgs,"Puccinia graminis")[[1]]

select(orgdb, head(keys(orgdb)),"SYMBOL")
head(keys(orgdb),"SYMBOL")
################################################################################

Mgenes <- read.table('276genes.txt', header = FALSE, stringsAsFactors = FALSE)$V1

Mgids <- mapIds(orgdb,
               keys = Mgenes,
               column = "GID",
               keytype = "SYMBOL",
               multiVals = "first")
Mgenes <- as.character(Mgenes)
Mgenes <- Mgenes[!is.na(Mgenes)]

ego <- enrichGO(Mgids,
  OrgDb = orgdb,
  pvalueCutoff  = 0.5
)

goplot(ego)

################################################################################
################################################################################
# PART2: Go back to the python script to edit the goatools output

################################################################################
# PART3: Make plot

library(ggplot2)
library(ggrepel)
library(stringr)
library(ggh4x)
library(viridis)

bubble <- function(d, title) {
  # Filter the data based on adjusted_p_value
  d <- d[d$adjusted_p_value < 1,]
  
  cols <- c("#4059AD", "#F4B942", "#b2df8a")
  
  bubble_plot <- ggplot(d, aes(x = term, y = -log10(adjusted_p_value), colour = source, size = intersection_size, label = term)) +
    geom_point(aes(fill = factor(source), alpha = 0.5), pch = 21) +
    scale_fill_manual(values = cols) + scale_color_manual(values=cols) +
    xlab("") + ylab("-log10 (Padj)") +
    theme_classic() +
    scale_size(range = c(1, 10)) +
    geom_text_repel(data = subset(d, -log10(adjusted_p_value) > 0), aes(label = stringr::str_wrap(term, 2)), size = 3) +
    ggtitle(title) + geom_hline(aes(yintercept = 1.5)) +
    facet_wrap(~ source, scales = "free_x") +
    theme(legend.position = "none", panel.grid.major = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),panel.spacing = unit(0, "lines"),
          strip.background = element_rect(color='white', fill='#97D8C4'))
  
  print(bubble_plot)
}

lollipop_plot <- function(d, title) {
  # Filter the data based on adjusted_p_value
  d <- d[d$adjusted_p_value < 1,]
  
  # Convert ratio_in_study to number
  d$ratio_in_study <- sapply(strsplit(d$ratio_in_study, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  
  lollipop_plot <- ggplot(d, aes(x = ratio_in_study, y = term, colour = -log10(adjusted_p_value), label = term)) +
    geom_segment(aes(yend = term, xend = 0), color = "grey25", alpha=0.5) +
    geom_point(aes(fill = -log10(adjusted_p_value), size = intersection_size), shape = 21, color = "grey25") +
    scale_fill_gradient("-log10(p-val)",low='#58A7A0', high='#7858A7') +
    scale_color_gradient("-log10(p-val)",low='#58A7A0', high='#7858A7') +
    ylab("") + xlab("Gene Ratio") +
    theme_minimal() + scale_x_continuous(limits=c(0,0.12), breaks=seq(0,0.12,0.03)) +
    scale_size("Count",range = c(1, 10)) +
    geom_text_repel(data = subset(d, ratio_in_study > 0), aes(label = stringr::str_wrap(term, 0)), size = 3) +
    ggtitle(title) + geom_vline(aes(xintercept=0)) +
    theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  print(lollipop_plot)
}

d <- read.csv("MARPLE_goenrich_Rinp.csv")
lollipop_plot(d,"")

ggsave("goenrich_plt_narrow.pdf",lollipop_plot(d,"GO-Enrichment of MARPLE genes"), width=4, height=5)

bubble(d, "Cluster # GO Enrichment Analysis")
