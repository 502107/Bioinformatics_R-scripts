library(adegenet)
library(reshape2)
library("adegenet")
library("stringr")
library(tibble)
library(scales)
library(ggplot2)

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

maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
set.seed(26)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(my_genind, n.pca = 75, max.n.clust = maxK, choose.n.clust = FALSE)
  myMat[i,] <- grp$Kstat
}

my_df2 <- melt(myMat)
colnames(my_df2)[1:3] <- c("Group", "K", "BIC")
my_df2$K <- as.factor(my_df2$K)
head(my_df2)

p1 <- ggplot(my_df2, aes(x = K, y = BIC)) +
  geom_boxplot(aes(fill = (K == 5)), lwd = 0.15, alpha=0.8) +
  scale_fill_manual(values = c("TRUE" = "palegreen1", "FALSE" = "white")) +
  theme_bw() +
  xlab("Number of clades (K)") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "cornflowerblue", linewidth=0.3, alpha=1) +
  geom_point(aes(color = (K == 5), shape=(K==5)), size=2, alpha=0.6) +
  scale_shape_manual(values=c("TRUE"=18,'FALSE'=20)) +
  scale_colour_manual(values = c("TRUE"='darkorchid1','FALSE'="black")) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
p1

#ggsave('all_BIC.pdf',p1, width = 16, height=8)
###############################################################################

my_k <- 13

grp_l <- vector(mode = 'list', length = length(my_k))
dapc_l <- vector(mode = 'list', length = length(my_k))

set.seed(1)
grp_l <- find.clusters(my_genind, n.pca=90, n.clust = 13)
dapc_l <- dapc(my_genind, pop = grp_l$grp, n.pca = 90, n.da = 4)

my_df2 <- as.data.frame(dapc_l$ind.coord)
my_df2$Group <- grp_l$grp
my_df2$index <- rownames(my_df2)

rename_isols <- c(
  "SAMN00988997" = "01MN84A-1-2",
  "SAMN00988993" = "06ND76C",
  "SAMN00988960" = "07KEN11-2",
  "SAMN00988968" = "07KEN24-4",
  "SAMN00988979" = "09ETH8-3",
  "SAMN00988959" = "09TAN06-2",
  "SAMN00988963" = "09TAN10-1",
  "SAMN00988967" = "09TAN8-16",
  "SAMN02887294" = "126-6711",
  "SRR9029856" = "194",
  "SRR6242031" = "21-0",
  "SRR6242032" = "326-1-2-3-5-6",
  "SRR6242048" = "34-2-12",
  "SRR6242033" = "34-2-12-13",
  "SAMN00989004" = "56SD37B",
  "SAMN00988992" = "59KS19",
  "SAMN00989000" = "59OH5B",
  "SAMN00988980" = "60IRN10B",
  "SAMN00989003" = "65SD2",
  "SAMN00989001" = "71MN603-2",
  "SAMN00988981" = "72ETH11-4",
  "SAMN00988996" = "74MN1409",
  "SAMN00988999" = "75MN68C",
  "SAMN00988994" = "77ND82A",
  "SAMN00988976" = "80MAR33B",
  "SAMN00988972" = "83ETH06-1",
  "SAMN00988974" = "83KEN7A",
  "SAMN00988982" = "84CSK759C",
  "SAMN00988990" = "84CSK764-3",
  "SAMN00988983" = "84ETH19A",
  "SAMN00988987" = "84KEN8C",
  "SAMN00988975" = "85MAR42-2",
  "SAMN00988978" = "86ITA1042A",
  "SAMN00988991" = "86MAR185C",
  "SAMN00988977" = "86PAK1030A",
  "SAMN00988973" = "87KEN11-4",
  "SAMN00988986" = "87KEN3018-4",
  "SAMN00988971" = "96ZIM2A",
  "SRR9029855" = "98",
  "SAMN00988995" = "99KS76A",
  "ERR2099117" = "CZ-01",
  "ERR2099118" = "CZ-02",
  "ERR2099119" = "CZ-03",
  "ERR2099114" = "DK-01",
  "ERR2099115" = "DK-02",
  "SAMN00988965" = "Eri-2010-13A",
  "ERR2099120" = "ET-01",
  "ERR2099121" = "ET-02",
  "ERR2099122" = "ET-03",
  "ERR2099123" = "ET-04",
  "ERR2099135" = "H-01",
  "ERR2099124" = "IR-01",
  "ERR2099125" = "IR-02",
  "ERR2099126" = "IR-03",
  "ERR2099127" = "IR-04",
  "ERR2099128" = "IR-05",
  "ERR2099129" = "IR-06",
  "ERR2099131" = "IS-02",
  "ERR2099132" = "IS-03",
  "SAMN00988984" = "ISR2147",
  "SAMN00988985" = "ISR2164",
  "ERR2099149" = "IT-01",
  "ERR2099150" = "IT-02",
  "ERR2099151" = "IT-03",
  "ERR2099152" = "IT-04",
  "SAMN00988966" = "KEN2010-10A",
  "ERR2099133" = "M-01",
  "SRR5883046" = "Pgt279",
  "SRR5883047" = "Pgt632",
  "ERR2099140" = "SA-01",
  "ERR2099141" = "SA-02",
  "ERR2099142" = "SA-03",
  "ERR2099143" = "SA-04",
  "ERR2099144" = "SA-05",
  "ERR2099145" = "SA-06",
  "ERR2099146" = "SA-07",
  "ERR2099116" = "SE-01",
  "ERR2099155" = "SE-02",
  "ERR2099156" = "SE-03",
  "SRR9024786" = "Ug99",
  "ERR2099153" = "UK-01",
  "ERR2099147" = "UR-01",
  "ERR2099148" = "UR-02",
  "ERR2099154" = "US-01",
  "SRR9024811" = "UVPgt55",
  "SRR9024812" = "UVPgt59",
  "SRR9024809" = "UVPgt60",
  "SRR9024810" = "UVPgt61",
  "SRR17888316" = "06KEN19V3",
  "SRR17888337" = "04KEN15604",
  "SRR17888351" = "84KEN8B",
  "SRR17888345" = "87KEN3018-1",
  "SRR17888344" = "87MDG1054A",
  "SRR17888350" = "84MAR10A",
  "SRR17888384" = "13GER17-2",
  "SRR17888387" = "13ETH22-2",
  "SRR17888379" = "69KS170",
  "SRR17888372" = "74MN1409-VI-A",
  "SRR17888393" = "00M063C",
  "SRR17888390" = "09ID073-2",
  "SRR17888364" = "76WA1433B",
  "SRR17888373" = "73WA399C",
  "SRR17888346" = "87ETH3008D",
  "SRR17888348" = "01TUR34A",
  "SRR17888381" = "01MN84-A-1-2",
  "SRR17888358" = "81BRA38A"
)

my_df2$index <- rename_isols[my_df2$index]

my_pal <- viridis::viridis(13)
my_pal2 <- viridis::viridis(13, direction=-1)

my_pal <- pals::cols25(13)

p2 <- ggplot(my_df2, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 5, shape = 21) +
  theme_bw() +
  scale_color_manual(values = my_pal) +
  scale_fill_manual(values = paste0(substr(my_pal, 1, 7), "66")) +
  labs(color = "Clade", fill = "Clade")
p2

p3 <- ggplot(my_df2, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_hex(bins = 1) +
  geom_jitter(size = 5, shape = 21, width = 0.3, height = 0.3) +
  theme_bw() +
  scale_color_manual(values = my_pal) +
  scale_fill_manual(values = paste0(substr(my_pal, 1, 7), "66")) +
  labs(color = "Clade", fill = "Clade")

p3





