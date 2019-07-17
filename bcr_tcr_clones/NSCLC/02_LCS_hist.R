library(SummarizedExperiment)
library(dplyr)
library(BioQC)
library(stringdist)
library(data.table)
library(BuenColors)
library(qualV)

# Filter down for consensus clonotypes
bcell <- read.table("../data/nsclc/vdj_v1_hs_nsclc_b_all_contig_annotations.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
tcell <- read.table("../data/nsclc/vdj_v1_hs_nsclc_t_all_contig_annotations.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

count_pair_LCS <- function(in_df, seed = 0){
  
  if(seed == 0){
    what <- "observed"
  } else{
    set.seed(seed)
    in_df$barcode <- sample(as.character(in_df$barcode))
    what <- "permuted"
  }
  
  in_df %>% filter(is_cell == "True"  & raw_clonotype_id != "None") %>%
    dplyr::select(barcode, raw_clonotype_id) %>% distinct() %>%
    group_by(raw_clonotype_id) %>%
    mutate(n = n()) %>% filter(n >= 2)  %>% arrange(raw_clonotype_id) -> df
  
  # Note: returns -inf if 0
  LCSCL <- function(a,b){
    
    aa <- strsplit(a, "")[[1]]
    bb <- strsplit(b, "")[[1]]
    
    rle_out <- rle(aa == bb)
    value <- max(rle_out$lengths[rle_out$values])
    value <- ifelse(value >0, value, 0)
    return(value)
  }
  
  # Determine Hamminng distance depending on degree of bead merging
  unique_clone_id <-df %>% pull(raw_clonotype_id) %>% unique()
  
  dist_df <- lapply(unique_clone_id, function(clone_name){
    
    # Determine appropriate barcodes
    barcodes <- gsub("-1", "", df %>% filter(raw_clonotype_id == clone_name) %>% pull(barcode) %>% as.character())
    if(length(barcodes) > 1){
      comb_mat <- combn(barcodes, 2)
      
      # Determine pair-wise Hamming distances
      dff <- lapply(1:dim(comb_mat)[2], function(i){
        dfo <- data.frame(barcode1 = comb_mat[1,i], barcode2 =  comb_mat[2,i], clone_name,
                          nbarcodes = length(barcodes),
                          LCS = LCSCL(comb_mat[1,i],  comb_mat[2,i]))
        dfo
      }) %>% rbindlist()  %>% data.frame()
    }
    dff 
  }) %>% rbindlist()
  
  # Summarize values and append zeros if necessary
  quant_df <- data.frame(value =  as.numeric(names(table(dist_df$LCS))), 
                         number = as.numeric(table(dist_df$LCS))/length(dist_df$LCS), what, seed)
  zero_df <- data.frame(
    value = 0:16,
    number = 0,
    what, seed
  )
  quant_df2 <- rbind(quant_df, zero_df) %>%
    group_by(what, seed, value) %>%
    summarize(number = sum(number))
  
  return(quant_df2)
}

dfb <- rbind(count_pair_LCS(bcell), count_pair_LCS(bcell, 1) )
dft <- rbind(count_pair_LCS(tcell), count_pair_LCS(tcell, 1) )

dfb %>% group_by(what) %>% summarize(sum(number))
dft %>% group_by(what) %>% summarize(sum(number))

dfb %>% filter(value >= 9) %>% group_by(what) %>%
  summarize(total = sum(number))
dft %>% filter(value >= 9) %>% group_by(what) %>%
  summarize(total = sum(number))

p1 <- ggplot(dfb, aes(x = value, y = number*100, fill = what)) +
  geom_bar(stat = "identity", position = "dodge", color = "black" , width = 0.8) +
  labs(x = "Restricted longest common subsequence", y = "% of barcode pairs within clonotype") +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = c(0.8, 0.8)) +
  scale_fill_manual(values = c("darkgrey", "lightgrey")) 


p2 <- ggplot(dfb %>% filter(value >=9), aes(x = value, y = number*100, fill = what)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.8) +
  labs(x = "Restricted longest common subsequence", y = "% of barcode pairs within clonotype") +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = c(0.8, 0.8)) +
  scale_fill_manual(values = c("darkgrey", "lightgrey")) 

cowplot::ggsave(cowplot::plot_grid(p1, p2, nrow = 1), file = "plots/rLCS_Bcell.pdf", width = 4, height = 2)

