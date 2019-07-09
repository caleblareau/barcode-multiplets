library(SummarizedExperiment)
library(dplyr)
library(BioQC)
library(stringdist)
library(data.table)
library(BuenColors)
library(qualV)

# Filter down for consensus clonotypes
tcell_clones <- read.table("../data/nsclc/vdj_v1_hs_nsclc_b_all_contig_annotations.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
tcell_clones %>% filter(is_cell == "True"  & raw_clonotype_id != "None") %>%
  dplyr::select(barcode, raw_clonotype_id) %>% distinct() %>%
  group_by(raw_clonotype_id) %>%
  mutate(n = n()) %>% filter(n >= 2)  %>% arrange(raw_clonotype_id) -> df

LCSCL <- function(a,b){
  rle_out <- rle(a == b)
  return(max(rle_out$lengths[rle_out$values]))
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
      data.frame(barcode1 = comb_mat[1,i], barcode2 =  comb_mat[2,i], clone_name,
                 nbarcodes = length(barcodes),
                 dist = stringdist(comb_mat[1,i], comb_mat[2,i]),
                 LCS = (LCSCL(strsplit(comb_mat[1,i], "")[[1]],  strsplit(comb_mat[2,i], "")[[1]])))
    }) %>% rbindlist()  %>% data.frame()
  }
  dff 
}) %>% rbindlist()

qplot(dist_df$LCS) +
  labs(x = "Longest common subsequence", y = "count") +
  pretty_plot() + L_border()
#dist_df %>% filter(clone_name == "clonotype11")
dist_df %>% filter(LCS >= 6)


tcell_clones %>% filter( raw_clonotype_id  == "clonotype9") %>%
  filter(is_cell == "True"  & raw_clonotype_id != "None") %>%
  pull("barcode") %>% unique() %>% data.frame() %>% write.table(quote = FALSE, row.names = FALSE)
