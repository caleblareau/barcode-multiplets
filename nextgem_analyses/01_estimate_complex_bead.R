library(dplyr)
library(data.table)
library(stringr)
library(BuenColors)
source("00_functions.R")

# Import barcode pairs
file <- "atac_pbmc_5k_v1_possorted_bam.barcodeTranslate.tsv"; idx <- 8
#file <- "atac_pbmc_5k_nextgem_possorted_bam.barcodeTranslate.tsv"; idx <- 8

dd <- fread(file, header = FALSE)
ss <- str_split_fixed(dd[[2]], "_", 8)
dd2 <- dd[ss[,idx] != "N01",] %>% arrange(V2)
dd2$barcode <- gsub("-1", "", dd2[,1])

# Get multiple rate
sum(ss[,8] != "N01")/dim(ss)[1]

# Loop through multiplets and compute the rLCS
rLCS <- function(a,b){
  
  aa <- strsplit(a, "")[[1]]
  bb <- strsplit(b, "")[[1]]
  
  rle_out <- rle(aa == bb)
  value <- max(rle_out$lengths[rle_out$values])
  value <- ifelse(value >0, value, 0)
  return(value)
}

rLCS_vec<- function(vec){
  mat <- combn(vec,2)
  sapply(1:dim(mat)[2], function(i){
    rLCS(mat[1,i], mat[2,i])
  }) %>% mean()
}

multiplets <- unique(dd2$V2)
means_vec <- sapply(multiplets, function(one_merged_bc){
  bcs <- dd2 %>% filter(V2 == one_merged_bc) %>% pull(barcode)
  ifelse(length(bcs) > 1, rLCS_vec(bcs), NA) # mostly a fail-safe; the NA shouldn't be thrown here...
})
multiplets_non_barcode_similarity <- multiplets[means_vec < 6 & !is.na(means_vec)]
multiplets_from_barcode_similarity <- multiplets[means_vec >= 6 & !is.na(means_vec)]

# Estimate the number of beads per droplet 
qdf <- data.frame(number = str_split_fixed(multiplets_non_barcode_similarity, "_", 8)[,idx]) %>%
  group_by(number) %>% summarize(count = n()) 
qdf2 <- data.frame(number = str_split_fixed(multiplets_from_barcode_similarity, "_", 8)[,idx]) %>%
  group_by(number) %>% summarize(count = n()) 

qdf$n <- as.numeric(gsub("N0", "", qdf$number)); qdf2$n <- as.numeric(gsub("N", "", qdf2$number))

p1 <- ggplot(qdf, aes(x = number, y = count)) +
  geom_bar(stat = "identity", color = "black", fill = "lightgrey") +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Number of beads / drop (from bap)", y = "count")

#cowplot::ggsave(p1, file = "../plots/R2R_beaddoublet_distribution.pdf", width = 2, height =2)

# Estimate rate of bead heterogenity
means_vec <- means_vec[!is.na(means_vec)]
n_heterogeneous_beads <- sum(means_vec >= 6)

n_barcodes_from_heterogeneous_beads <- sum(qdf2$count)

# double count physical doublets but remove extras from heterogeneous beads
n_total_beads <- dim(dd)[1] - sum((qdf2$n -1) * (qdf2$count))
n_heterogeneous_beads/n_total_beads

# propotrion that are from complex bead
sum(means_vec >= 6)/length(means_vec)
