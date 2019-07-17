library(dplyr)
library(BuenColors)
library(stringr)
library(stringdist)
library(data.table)
library(ggrastr)

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

set.seed(6)
idx <- sample(1:dim(tenx_5k)[1])
p0 <- ggplot(tenx_5k[idx,], aes(x = TSNE.1, y = TSNE.2, color = as.factor(Cluster))) +
  geom_point_rast(dpi = 1000, size = 0.1) + 
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "t-SNE 1", y = "t-SNE 2" ) 
cowplot::ggsave(p0, file = "../plots/tSNE_basic.pdf", width = 1.8, height = 1.8)


set.seed(6)
idx <- sample(1:dim(tenx_5k)[1])
p1 <- ggplot(tenx_5k[idx,], aes(x = TSNE.1, y = TSNE.2, color = isMultiplet)) +
  geom_point_rast(dpi = 500, size = 0.1, alpha = 0.5) + 
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("lightgrey", "red"))
cowplot::ggsave(p1, file = "../plots/tSNE_5k_out.pdf", width = 2, height = 2)

zoom_df <- tenx_5k %>% filter(TSNE.1 > 25 & TSNE.2 < -10)
zoom_df_plot <- zoom_df %>% group_by(bap_id) %>% mutate(n = n()) %>%
  mutate(bid = as.character(bap_id)) %>%
  ungroup() %>%
  mutate(color = ifelse(n == 1, "Singlet", as.character(bap_id)))


p2 <- ggplot(zoom_df_plot %>% arrange(desc(color)), aes(x = TSNE.1, y = TSNE.2, color = color)) +
  geom_point(size = 1) + 
  pretty_plot() + L_border() +
  scale_color_manual(values = c("purple2", "firebrick", "green3", "dodgerblue3", "orange3", "grey")) +
  theme(legend.position = "right") +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank())

cowplot::ggsave(p2, file = "../plots/tSNE_zoom.pdf", width = 1.5, height = 1.8)




p3 <- ggplot(tenx_5k %>% arrange(desc(bc9)), aes(x = TSNE.1, y = TSNE.2, color = bc9)) +
  geom_point_rast(dpi = 1000, size = 0.1) + 
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") + labs(x = "t-SNE1", y = "t-SNE2" ) +
  scale_color_manual(values = c("firebrick", "dodgerblue", "lightgrey"))

cowplot::ggsave(p3, file = "../plots/tSNE_bc9.pdf", width = 1.8, height = 1.8)


