library("adegenet")
library("stringr")
library(tibble)
library(ggplot2)
library(reshape2)
library(pals)
library(dplyr)

df = read.table("structure_86isols.txt", colClasses = "character")
individuals = unique(str_sub(colnames(df), 1, -3))

# Convert from 2 column per locus to 1 column per locus
for (individual in individuals) {
  allele_one_col = paste(individual, "_1", sep="")
  allele_two_col = paste(individual,"_2", sep="")
  df[[individual]] = paste(df[[allele_one_col]], df[[allele_two_col]])
}

# Keep only the columns that hold both alleles
df = df[individuals]
df = t(df)
my_genind = df2genind(df, sep = " ", NA.char = "-9", check.ploidy = F)

# Free up memory by removing original dataframe
rm(df)

# Get number of loci used for genind object
nLoc(my_genind)

set.seed(3)
my_clus = find.clusters(my_genind, n.pca = 90, max.n.clust = 15)
my_dapc = dapc(my_genind, my_clus$grp, n.pca=90, n.da=4)
my_dapc

write.csv(cbind(my_dapc$assign, my_dapc$posterior), file="assignations104_spltDAPC/assign_90_7_90_4_seed8.csv")

################################################################################
# Basic scatter plot with isolates
myCol <- c("darkblue","purple","green","orange","red","blue","pink","darkgreen","magenta")
scatter(my_dapc, col=myCol, cex=1.5, clabel=0, posi.da="bottomright")
scatter(my_dapc,1,2, col=my_dapc$grp, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)
compoplot(my_dapc, subset=1:50, posi="bottomright",
          txt.leg=paste("Cluster", 1:7), lab="",
          ncol=2, xlab="individuals", col=funky(6))
scatter(my_dapc, grp = my_dapc$grp)
scatter(my_dapc, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")

cols <- c(
  '#0096c1',
  '#877192',
  '#fac8c8',
  '#e5a0bf',
  '#fff1bc',
  '#e89b5c',
  '#add3bb'
  )

scatter(my_dapc, scree.da=FALSE, bg="white", pch=20, cell=1, cstar=1, col=cols, solid=1, cex=0.1,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:7))
scatter.dapc(my_dapc,scree.pca=FALSE, scree.da=FALSE, legend=TRUE)

################################################################################
# Extract eigenvalues from PCA step
eigenvalues <- my_dapc$pca.eig

if(length(eigenvalues) > 0) {
  variation_percentage <- (eigenvalues / sum(eigenvalues)) * 100
  variation_df <- data.frame(PC = 1:length(variation_percentage), 
                             Variation = variation_percentage)
  
  # Plot the percentage of variation
  ggplot(variation_df, aes(x = PC, y = Variation)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(x = "PC",
         y = "Variation (%)") +
    geom_text(aes(label = round(Variation, 2)), vjust = -0.5, size = 3)
} else {
  print("No eigenvalues found in PCA results.")
}

da_eigenvalues <- as.numeric(my_dapc$eig)
if(length(da_eigenvalues) > 0) {
  # Calculate the percentage of variation for each discriminant function
  da_variation_percentage <- (da_eigenvalues / sum(da_eigenvalues)) * 100
  da_variation_df <- data.frame(DiscriminantFunction = 1:length(da_variation_percentage), 
                                Variation = da_variation_percentage)
  
  # Plot the percentage of variation
  ggplot(da_variation_df, aes(x = DiscriminantFunction, y = Variation)) +
    geom_bar(stat = "identity", fill = "coral") +
    theme_minimal() +
    labs(title = "Percentage of Variation Captured by Each Discriminant Function",
         x = "Discriminant Function",
         y = "Percentage of Variation (%)") +
    geom_text(aes(label = round(Variation, 2)), vjust = -0.5, size = 3)
} else {
  print("No eigenvalues found in discriminant analysis results.")
}

# Extract posterior probabilities
posterior_probs <- as.data.frame(my_dapc$posterior)
posterior_probs$Individual <- rownames(posterior_probs)
posterior_probs_melt <- melt(posterior_probs, id.vars = "Individual", variable.name = "Cluster", value.name = "Probability")

# Plot posterior membership probabilities
ggplot(posterior_probs_melt, aes(x = Cluster, y = Probability, fill = Cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Posterior Membership Probabilities for Each Cluster",
       x = "Cluster",
       y = "Posterior Probability")

# Extract coordinates of individuals in PCA space
custom_cl <- c(1,4,5,7,2,3,6)
pca_coords <- as.data.frame(my_dapc$ind.coord)
pca_coords$Cluster <- my_dapc$assign
cluster_variation <- pca_coords %>%
  group_by(Cluster) %>%
  summarize(across(everything(), var))

cluster_variation_melt <- melt(cluster_variation, id.vars = "Cluster", variable.name = "PC", value.name = "Variance")

# Plot within-cluster variation
ggplot(cluster_variation_melt, aes(x = PC, y = Variance, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Variation for Each Principal Component",
       x = "PC",
       y = "Variance") +
  scale_fill_manual(values = cols)


################################################################################

assignplot(my_dapc)

df <- cbind(my_dapc$assign,my_dapc$posterior)
df <- as.data.frame(df)

df <- rownames_to_column(df, var = "accession")
names(df)[names(df) == 'V1'] <- "cluster"

order_names <- read.table(sep='\t', file='edits/acc_treenames_clades.tsv')
colnames(order_names) <- c('accession','names','order','clade')
df_merge <- merge(df,order_names,by="accession")

df_long <- melt(df_merge, id.vars = c("accession","names","cluster","order","clade"))

ggplot(df_long, aes(x = order, y = value, fill = variable)) +
  geom_bar(stat = "identity", width=1) +
  theme_minimal() + scale_y_continuous(expand=c(0,0)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  labs(x = "Accession", y = "Percentage", fill = "Clade") +
  scale_fill_manual(values=as.vector(stepped(20)))

################################################################################
# Add the clade plot next to the phylogenetic tree

df_txt <- read.table("edits/new_Mclades_accessionnames.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df_long <- merge(df_long, df_txt, by = "accession", all.x = TRUE)
df_long$accession <- factor(df_long$accession, levels = df_long$accession[order(df_long$order)])
write.csv(df_long,"check.csv")

p1 <- ggplot(df_long, aes(x = order, y = value, fill = variable)) +
  geom_bar(stat = "identity", width=1) +
  theme_minimal() + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  labs(x = "Accession", y = "Percentage", fill = "Clade") +
  scale_fill_manual(values=as.vector(stepped(20)))
  

ggsave("ordered_clades_2.pdf",plot=p1)
dev.off()
