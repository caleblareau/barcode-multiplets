library(dplyr)
library(data.table)
library(stringr)
library(BuenColors)

# Import barcode pairs
file <- "../data/pub_5k/atac_v1_pbmc_5k_possorted_bam.barcodeTranslate.tsv"; idx <- 8
file <- "../data/sai/possorted_bam.barcodeTranslate.tsv"; idx <- 4

dd <- fread(file, header = FALSE)
ss <- str_split_fixed(dd[[2]], "_", 8)
dd2 <- dd[ss[,idx] != "N01",] %>% arrange(V2)
dd2$barcode <- gsub("-1", "", dd2[,1])

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

qdf$n <- as.numeric(gsub("N0", "", qdf$number)); qdf2$n <- as.numeric(gsub("N0", "", qdf2$number))

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
n_total_beads <- dim(dd)[1] + sum((qdf$n -1) * (qdf$count)) - sum((qdf2$n -1) * (qdf2$count))
n_heterogeneous_beads/n_total_beads

sum(means_vec >= 6)/length(means_vec)

# Compare the abundances between the types of multiplets
pct_diff <- function(a,b){
  return(abs(a-b)/mean(a,b)*100)
}

pct_vec <-  function(vec){
  mat <- combn(vec,2)
  sapply(1:dim(mat)[2], function(i){
    pct_diff(mat[1,i], mat[2,i])
  }) %>% mean()
}

sc <- fread("../data/pub_5k/atac_v1_pbmc_5k_singlecell.csv") %>% filter(cell_id != "None") %>% data.f

compute_percent_differences_vector <- function(bc_groups){
  vec_out <- sapply(bc_groups, function(bc_group){
    bcs <- dd %>% filter(V2 == bc_group) %>% pull(V1) %>% as.character
    vec <- sc %>% filter(barcode %in% bcs) %>% pull(passed_filters) %>% log2
    pct_vec(vec)
  })
  return(vec_out)
}
multibarcode_1bead <- compute_percent_differences_vector(multiplets_from_barcode_similarity)
multibead <- compute_percent_differences_vector(multiplets_non_barcode_similarity)
pct_df <- data.frame(
  value = c(multibarcode_1bead, multibead), 
  type = c(rep("B", length(multibarcode_1bead)), rep("A", length(multibead)))
)
ks.test(multibarcode_1bead, multibead)
p00 <- ggplot(pct_df, aes(x = type, y = value, color = type)) + 
  geom_boxplot(width = 0.3, fill = NA, outlier.shape = NA)+
  scale_color_manual(values = c("black", "firebrick")) +
   labs(x = "Inferred multiplet type", y = "Mean % difference of log2 fragments") +
  pretty_plot(fontsize = 8) + L_border()+
  theme(legend.position = "none") 
cowplot::ggsave(p00, file = "../plots/pct_diff_amount.pdf", width = 1.3, height = 2)

# Plot the distance
p0 <- ggplot(data.frame(means_vec), aes(x = means_vec)) +
  geom_density(adjust = 4, fill = "lightgrey") +
  labs(x = "Mean rLCS per multiplet", y = "density") +
  pretty_plot(fontsize = 8) + scale_y_continuous(expand = c(0,0)) +
  L_border() + geom_vline(xintercept = 5.8, linetype = 2, color = "firebrick")

#cowplot::ggsave(p0, file = "../plots/rLCS_density.pdf", width = 2, height = 2)
