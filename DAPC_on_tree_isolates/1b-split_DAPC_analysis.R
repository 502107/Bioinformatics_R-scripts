library(adegenet)
library(stringr)
library(tibble)
library(ggplot2)
library(reshape2)
library(pals)
library(dplyr)
library(ape)
library(ggtree)

###############################################################################

df = read.table("structure_86isols.txt", colClasses = "character")
individuals = unique(str_sub(colnames(df), 1, -3))
for (individual in individuals) {
  allele_one_col = paste(individual, "_1", sep="")
  allele_two_col = paste(individual,"_2", sep="")
  df[[individual]] = paste(df[[allele_one_col]], df[[allele_two_col]])
}

df = df[individuals]
df = t(df)
my_genind = df2genind(df, sep = " ", NA.char = "-9", check.ploidy = F)
rm(df)
nLoc(my_genind)

###############################################################################

set.seed(10)
my_clus = find.clusters(my_genind, n.pca = 75, max.n.clust = 15)

###############################################################################
my_dapc = dapc(my_genind, my_clus$grp, n.pca=75, n.da=4)
#write.csv(cbind(my_dapc$assign, my_dapc$posterior), file="assignation_results/assign_90_15_90_4_seed10.csv")

df <- cbind(my_dapc$assign,my_dapc$posterior)
df <- as.data.frame(df)

df <- rownames_to_column(df, var = "accession")
names(df)[names(df) == 'V1'] <- "cluster"

order_names <- read.table(sep='\t', file='edits/old_acc_treenames_clades.tsv')
colnames(order_names) <- c('accession','names','order','clade')
df_merge <- merge(df,order_names,by="accession")

df_long <- melt(df_merge, id.vars = c("accession","names","cluster","order","clade"))

###############################################################################
# Plot DAPC clusters

cols <- c(
  '#6baed6', 
  '#9ecae1', 
  '#e6550d', 
  '#e6a8a8', 
  '#fd8d3c', 
  '#e8d5a1', 
  '#c6dbef', 
  '#31a354', 
  '#74c476', 
  '#a1d99b', 
  '#c7e9c0', 
  '#c7e9c0', 
  '#c7e9c0', 
  '#c87a41', 
  '#a4c0ce'  
)

cols5 <- c(
  '#88ccee',
  '#e57564',
  '#44aa99',
  '#e8d5a1',
  '#e6a8a8' 
)


ggplot(df_long, aes(x = order, y = value, fill = variable)) +
  geom_bar(stat = "identity", width=1) +
  theme_minimal() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  labs(x = "Accession", y = "Percentage", fill = "Clade") +
  scale_fill_manual(values=as.vector(cols5))

###############################################################################
# Select sub-clade to run DAPC again

df2 = read.table("structure_subclade_inp/Structure1.txt", colClasses = "character")
individuals = unique(str_sub(colnames(df2), 1, -3))
for (individual in individuals) {
  allele_one_col = paste(individual, "_1", sep="")
  allele_two_col = paste(individual,"_2", sep="")
  df2[[individual]] = paste(df2[[allele_one_col]], df2[[allele_two_col]])
}

df2 = df2[individuals]
df2 = t(df2)
my_genind2 = df2genind(df2, sep = " ", NA.char = "-9", check.ploidy = F)
rm(df2)
nLoc(my_genind2)

set.seed(8)
my_clus2 = find.clusters(my_genind2, n.pca = 45, max.n.clust = 10)
my_dapc2 = dapc(my_genind2, my_clus2$grp, n.pca=4, n.da=4)
#write.csv(cbind(my_dapc2$assign, my_dapc2$posterior), file="subclade_assignation_results/assign_pc45_cl3_sd26_oo10.csv")

df2 <- cbind(my_dapc2$assign,my_dapc2$posterior)
df2 <- as.data.frame(df2)

df2 <- rownames_to_column(df2, var = "accession")
names(df2)[names(df2) == 'V1'] <- "cluster"

df_merge2 <- merge(df2,order_names,by="accession")
df_long2 <- melt(df_merge2, id.vars = c("accession","names","cluster","order","clade"))

p <- ggplot(df_long2, aes(x = order, y = value, fill = variable)) +
  geom_bar(stat = "identity", width=1) +
  theme_minimal() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  labs(x = "Accession", y = "Percentage", fill = "Clade") +
  scale_fill_manual(values=tail(stepped3(20),-4))

p



#ggsave('Stucture2_3cl_stepped3-4.pdf',p)

                    