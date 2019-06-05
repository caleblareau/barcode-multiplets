library(dplyr)
library(BuenColors)
library(stringr)
library(stringdist)
library(data.table)

# Sai's data
cdf <- read.table("../data/sai/clusters.csv", sep = ",", header = TRUE)
bdf <- read.table("../data/sai/possorted_bam.barcodeTranslate.tsv"); colnames(bdf) <- c("Barcode", "bap_id")
bdf$n_beads <- str_sub(bdf$bap_id,-2,-1) %>% as.character() %>% as.numeric
tsnedf <- read.table("../data/sai/projection.csv", sep = ",", header = TRUE)
sai_df <- merge(merge(bdf, cdf), tsnedf)
sai_df$isMultiplet <- sai_df$n_beads > 1

ggplot(sai_df %>% arrange((isMultiplet)), aes(x = TSNE.1, y = TSNE.2, color = isMultiplet)) +
  geom_point() + 
  pretty_plot() + L_border() +
  scale_color_manual(values = c("black", "red"))


#--------------------------------------------

# 10X 5k data
cdf <- read.table("../data/pub_5k/clusters.csv", sep = ",", header = TRUE)
bdf <- read.table("../data/pub_5k/atac_v1_pbmc_5k_possorted_bam.barcodeTranslate.tsv"); colnames(bdf) <- c("Barcode", "bap_id")
bdf$n_beads <- str_sub(bdf$bap_id,-2,-1) %>% as.character() %>% as.numeric
tsnedf <- read.table("../data/pub_5k/projection.csv", sep = ",", header = TRUE)
tenx_5k <- merge(merge(bdf, cdf), tsnedf)
tenx_5k$isMultiplet <- tenx_5k$n_beads > 1
tenx_5k$bc9 <- ifelse(tenx_5k$bap_id == "atac_v1_pbmc_5k_possorted_bam_BC3848_N09", "group1", ifelse(tenx_5k$bap_id == "atac_v1_pbmc_5k_possorted_bam_BC2220_N09", "group2", "other"))

ggplot(tenx_5k %>% arrange((isMultiplet)), aes(x = TSNE.1, y = TSNE.2, color = isMultiplet)) +
  geom_point() + 
  pretty_plot() + L_border() +
  scale_color_manual(values = c("black", "red"))

ggplot(tenx_5k %>% arrange((isMultiplet)), aes(x = TSNE.1, y = TSNE.2, color = bc9)) +
  geom_point() + 
  pretty_plot() + L_border() +
  scale_color_manual(values = c("red", "dodgerblue", "grey"))

tenx_5k %>% filter(TSNE.1 > 25 & TSNE.2 < -10) %>%
  arrange(bap_id)

tenx_5k %>% filter(n_beads > 8)
