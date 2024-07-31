library(ggplot2)
library(reshape2)
library(pals)
library(dplyr)
library(tidyverse)
library(stringr)

files <- list.files(path = "assignation_results", pattern = "assign_90.*_90_.*\\.csv$", full.names = TRUE)

my_dapc <- map(files, ~ {
  df <- read.csv(.x, header = TRUE)
  colnames(df)[1:2] <- c("accession","V1")
  file_num <- str_extract(.x, "(?<=_90_)\\d+(?=_90_)")
  df$file <- paste0("k=", file_num)
  df$filen <- as.integer(file_num)
  return(df)
})

df_long <- lapply(my_dapc, function(df) {
  df <- melt(df, id.vars = c("V1","accession","file", "filen"))
  df$variable <- as.integer(str_replace(df$variable, "X", ""))
  return(df)
})

df_combined <- do.call(rbind, df_long)
df_combined$variable <- as.factor(df_combined$variable)
df_combined$file_factor <- factor(df_combined$file, levels = unique(df_combined$file[order(df_combined$filen)]))

df_combined <- df_combined |>
  mutate(file_accession = interaction(file_factor, accession, lex.order = TRUE),
         file_accession = fct_reorder(file_accession, V1))

order_names <- read.table(sep='\t', file='edits/old_acc_treenames_clades.tsv')
colnames(order_names) <- c('accession','names','order','clade')
df_merge <- merge(df_combined,order_names,by="accession")

cols <- c(
  '#0096c1',  # Original
  '#e5a0bf',  # Original
  '#fff1bc',  # Original
  '#add3bb',  # Original
  '#877192',  # Original
  '#fac8c8',  # Original
  '#e89b5c',  # Original
  '#007ea4',  # Darker shade of #0096c1
  '#d089a3',  # Darker shade of #e5a0bf
  '#e8d5a1',  # Darker shade of #fff1bc
  '#8fb3a0',  # Darker shade of #add3bb
  '#6d5a74',  # Darker shade of #877192
  '#e6a8a8',  # Darker shade of #fac8c8
  '#c87a41',  # Darker shade of #e89b5c
  '#a4c0ce'   # Complementary blend of the palette
)

cols20 <- c(
  '#88ccee', #1
  '#e57564', #2
  '#44aa99', #3
  '#e8d5a1', #4
  '#e6a8a8', #5
  '#332288', #6
  '#117733', #7
  '#999933', #8
  '#ddcc77', #9
  '#cc6677', #10
  '#882255', #11
  '#aa4499', #12
  '#44bbbb', #13
  '#eedd88', #14
  '#ffaabb', #15
  '#ee8866', #16
  '#56b4e9', #17
  '#009e73', #18
  '#f0e442', #19
  '#0072b2'  #20
)

plot <- ggplot(df_merge, aes(x = order, y = value, fill = variable)) +
  geom_bar(stat = "identity", width=1) +
  theme_minimal() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        panel.spacing = unit(0.05, "lines"), text=element_text(size=18)) +
  labs(x = "Accession", y = "Percentage", fill = "Cluster") +
  scale_fill_manual(values=as.vector(cols20)) +
  facet_wrap(~ file_factor, ncol = 1, strip.position = "right", scales = "free", dir = "v")

plot

ggsave("clusters_of_all_clades_seed10_2.pdf",plot, width=16,height=8)
