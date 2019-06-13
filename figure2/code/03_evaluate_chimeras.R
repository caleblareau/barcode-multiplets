options(stringsAsFactors = FALSE)
library(dplyr)
library(BuenColors)
library(stringr)
library(stringdist)
library(data.table)

# Sai's data
bdf <- read.table("../data/pub_5k/atac_v1_pbmc_5k_possorted_bam.barcodeTranslate.tsv"); colnames(bdf) <- c("Barcode", "bap_id")
barcodes <- bdf %>% pull(Barcode) %>% as.character()
barcodes <- gsub("-1", "", barcodes)
barcode_trans_vec <-  bdf %>% pull(bap_id) %>% as.character(); names(barcode_trans_vec) <- barcodes

# Import chimeras
chimera_df <- read.table("../data/pbmc_5k_atac_chimera.tsv"); colnames(chimera_df) <- c("parsed_barcode", "observed_barcode")
chimera_df_filt <- chimera_df %>% filter(parsed_barcode %in% barcodes & observed_barcode %in% barcodes)
chimera_df_filt_diff <- chimera_df_filt %>% filter(parsed_barcode != observed_barcode) 

# Annotate with bap barcodes
chimera_df_filt_diff$bap_1 <- barcode_trans_vec[as.character(chimera_df_filt_diff$parsed_barcode)]
chimera_df_filt_diff$bap_2 <- barcode_trans_vec[as.character(chimera_df_filt_diff$observed_barcode)]
chimera_df_filt_diff
