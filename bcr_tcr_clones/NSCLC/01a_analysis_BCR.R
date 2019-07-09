library(SummarizedExperiment)
library(dplyr)
library(BioQC)

# Filter down for consensus clonotypes
bcell_clones <- read.table("../data/nsclc/vdj_v1_hs_nsclc_b_all_contig_annotations.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
bcell_clones %>% filter(is_cell == "True"  & raw_clonotype_id != "None") %>%
  dplyr::select(barcode, raw_clonotype_id) %>% distinct() %>%
  group_by(raw_clonotype_id) %>%
  mutate(n = n()) %>% filter(n >= 3)  %>% arrange(raw_clonotype_id) -> filtered_b

# Identify HQ mutations
SE <- readRDS("../../NSCLC_mgatk/final/NSCLC-mgatk.rds")

computeAFMutMatrix <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), ref_allele, ">", letter)
    return(mat[ref_allele != letter,])
  }
  
  rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  
}

i <- intersect(as.character(filtered_b$barcode), colnames(SE))
af2 <- data.matrix(computeAFMutMatrix(SE[,i]))
clones <- filtered_b %>% filter(barcode %in% i) %>% pull(raw_clonotype_id)
af2 <- af2[rowMeans(af2) > 0.0001,]

lapply(names(table(clones)[table(clones) > 1]), function(clone){
  message(clone)
  # Do a Mann-Whitney test
  out <- sort(wmwTest(data.matrix(t(af2)), which(clones == clone),
                      valType = "abs.log10p.two.sided"), decreasing = TRUE)
  
  # Retain if reasonable p-value-- will threshold harder, later
  keep <- out[out> 4]
  c_1p <- rowSums(af2[names(keep),clones == clone] >= 0.01)
  nc_1p <- rowSums(af2[names(keep),clones != clone] >= 0.01)
  
  df <- data.frame(Variant = names(keep), log10p = unname(keep), colony = clone,
                   n_clones = sum(clones==clone), n_clone_1p = c_1p,
                   n_other_1p = nc_1p)
  df
})  %>% rbindlist() %>% as.data.frame() -> variantScoresClones

variantScoresClones %>% filter(n_clone_1p >= 3)


