library(data.table)
library(dplyr)
library(stringr)

# Import the multiplets calls from both sources; from the 10x excluded_barcodes, barcode 1 is what they filter
tenx_v1 <- fread("data/public/excluded_barcodes.csv", col.names = c("barcode1", "barcode2", "reason"))
bap_v1 <- fread("data/public/atac_v1_pbmc_5k_possorted_bam.barcodeTranslate.tsv", col.names = c("barcode", "bap_barcode"), header = FALSE)
passknee_barcodes_v1 <- bap_v1[["barcode"]]

# 10x reports all possible multiplets, not jsut those that pass the HQ barcodes, which is what bap considers
# So, filter the 10x data frame down such that the multiplets include the HQ barcodoes
tenx_v1_filt <- tenx_v1 %>% filter(barcode1 %in% passknee_barcodes_v1 & barcode2 %in% passknee_barcodes_v1)

# Now, we can define the implicated barcodes from 10x's solution
tenx_v1_implicated <- unique(c(tenx_v1_filt[["barcode1"]], tenx_v1_filt[["barcode2"]]))

# In bap, we define a barcode multiplet when the resulting merged multiplet has > 1 barcode
n_barcodes_multiplet <- as.numeric(gsub('N', '', str_split_fixed(bap_v1[["bap_barcode"]], "_", 8)[,8]))

# Now let's see per barcode how the two tools compare
df <- data.frame(
  barcode = passknee_barcodes_v1,
  bap_barcode = bap_v1[["bap_barcode"]],
  multiplet_bap = n_barcodes_multiplet > 1,
  multiplet_10x = passknee_barcodes_v1 %in% tenx_v1_implicated,
  n_barcodes_multiplet
)

#Rates
sum(df$multiplet_bap)/length(df$multiplet_bap)
sum(df$multiplet_10x)/length(df$multiplet_10x)

table(df$multiplet_bap, df$multiplet_10x)

# Concordance
sum(df$multiplet_bap == df$multiplet_10x)/dim(df)[1]

# These are likely 10x false negatives
df %>% filter(bap_barcode == "atac_v1_pbmc_5k_possorted_bam_BC4341_N04")
df %>% filter(bap_barcode == "atac_v1_pbmc_5k_possorted_bam_BC1157_N06")

# Here's probably a pair that 10x completely missed
df %>% filter(bap_barcode == "atac_v1_pbmc_5k_possorted_bam_BC0698_N02")

# Maybe this is a bap false positive
df %>% filter(bap_barcode == "atac_v1_pbmc_5k_possorted_bam_BC3438_N03")


# Import the multiplets calls from both sources; from the 10x excluded_barcodes, barcode 1 is what they filter
tenx_v1 <- fread("data/this_study/excluded_barcodes.csv", col.names = c("barcode1", "barcode2", "reason"))
bap_v1 <- fread("data/this_study/possorted_bam.barcodeTranslate.tsv", col.names = c("barcode", "bap_barcode"), header = FALSE)
passknee_barcodes_v1 <- bap_v1[["barcode"]]

# 10x reports all possible multiplets, not jsut those that pass the HQ barcodes, which is what bap considers
# So, filter the 10x data frame down such that the multiplets include the HQ barcodoes
tenx_v1_filt <- tenx_v1 %>% filter(barcode1 %in% passknee_barcodes_v1 & barcode2 %in% passknee_barcodes_v1)

# Now, we can define the implicated barcodes from 10x's solution
tenx_v1_implicated <- unique(c(tenx_v1_filt[["barcode1"]], tenx_v1_filt[["barcode2"]]))

# In bap, we define a barcode multiplet when the resulting merged multiplet has > 1 barcode
n_barcodes_multiplet <- as.numeric(gsub('N', '', str_split_fixed(bap_v1[["bap_barcode"]], "_", 4)[,4]))

# Now let's see per barcode how the two tools compare
df <- data.frame(
  barcode = passknee_barcodes_v1,
  bap_barcode = bap_v1[["bap_barcode"]],
  multiplet_bap = n_barcodes_multiplet > 1,
  multiplet_10x = passknee_barcodes_v1 %in% tenx_v1_implicated,
  n_barcodes_multiplet
)

# Rates
sum(df$multiplet_bap)/length(df$multiplet_bap)
sum(df$multiplet_10x)/length(df$multiplet_10x)

table(df$multiplet_bap, df$multiplet_10x)

# Concordance
sum(df$multiplet_bap == df$multiplet_10x)/dim(df)[1]


library(dplyr)
library(BuenColors)

df_plot <- data.frame(
  source = c("Public", "Public", "ThisStudy", "ThisStudy"), 
  what = c("bap", "10x", "bap", "10x"), 
  rate = c(17.6, 17.8, 13.2, 13.5)
)

p1 <- ggplot(df_plot, aes(x = source, y = rate, fill = what)) +
  geom_bar(color = "black", stat = "identity", position = position_dodge2(), width = 0.8) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values = c("darkgrey", "lightgrey"))+
  labs(x = "", y = "% barcode multiplets") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 24))  +
  theme(legend.position = "none")

cowplot::ggsave(p1, file = "rates_comparison.pdf", width = 2, height = 2)

