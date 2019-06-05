library(dplyr)
library(BuenColors)
library(stringr)
library(stringdist)
library(data.table)


# Determine number of bead barcodes that belong to a multiplet
prop_bead_barcodes_in_multiplet <- function(pub_df){
  table(pub_df$n_beads)/length((pub_df$n_beads))
}


# Determine number of droplets with multiplets
prop_droplets_multiplet <- function(pub_df){
  summary_df <- pub_df %>% group_by(bap_id) %>%
    summarize(n_beads_one = max(n_beads))
  table(summary_df$n_beads_one)/length(summary_df$n_beads_one)
}

# Chi-square tests
chisquare_cluster <- function(pub_df){
  pub_chisq_mat <- pub_df %>% mutate(Cluster = as.factor(as.character(Cluster)),
                                     multiplet = n_beads > 1) %>% 
    group_by(multiplet, Cluster) %>%
    summarize(count = n()) %>%
    reshape2::dcast(Cluster~multiplet, value.var = "count", fill = 0) %>%
    dplyr::select(c("FALSE", "TRUE")) %>% data.matrix() 
  
  chisq.test(pub_chisq_mat)
}

# Determine Hamminng distance depending on degree of bead merging
hamming_distance_count <- function(pub_df, what){
  unique_bap_pub <- pub_df %>% filter(n_beads > 1) %>% pull(bap_id) %>% as.character %>% unique
  nbead_dist_df <- lapply(unique_bap_pub, function(bap_name){
    
    # Determine appropriate barcodes
    barcodes <- pub_df %>% filter(bap_id == bap_name) %>% pull(Barcode) %>% as.character()
    if(length(barcodes) > 1){
      comb_mat <- combn(barcodes, 2)
      nbeads <- str_sub(bap_name,-2,-1)
      
      # Determine pair-wise Hamming distances
      dist_vec <- sapply(1:dim(comb_mat)[2], function(i){
        stringdist(comb_mat[1,i], comb_mat[2,i])
      })
      data.frame(
        dist_vec, 
        n_barcodes = as.character(as.numeric(as.character(nbeads))), stringsAsFactors = FALSE
      )
    }
    
  }) %>% rbindlist()
  
  nbead_dist_df %>% 
    group_by(n_barcodes) %>%
    summarize(mean = mean(dist_vec),
              sem = sqrt(var(dist_vec))/sqrt(n())) %>%
    mutate(what = what)
}

# Shared cluster for multiplet
prop_cluster_shared_multiplet <- function(pub_df, permute.seed = 0){
  
  if(permute.seed == 0){
    cluster_vec <- as.character(pub_df$Cluster); names(cluster_vec) <- as.character(pub_df$Barcode)
    
  } else{
    set.seed(permute.seed)
    cluster_vec <- sample(as.character(pub_df$Cluster)); names(cluster_vec) <- as.character(pub_df$Barcode)
    
  }
  unique_bap_pub <- pub_df %>% filter(n_beads > 1) %>% pull(bap_id) %>% as.character %>% unique
  same_cluster_df <- lapply(unique_bap_pub, function(bap_name){
    
    # Determine appropriate barcodes
    barcodes <- pub_df %>% filter(bap_id == bap_name) %>% pull(Barcode) %>% as.character()
    if(length(barcodes) > 1){
      comb_mat <- combn(barcodes, 2)
      nbeads <- str_sub(bap_name,-2,-1)
      
      # Determine pair-wise Hamming distances
      same_cluster <- sapply(1:dim(comb_mat)[2], function(i){
        as.numeric(cluster_vec[comb_mat[1,i]] == cluster_vec[comb_mat[2,i]])
      })
      data.frame(
        same_cluster, 
        n_barcodes = as.character(as.numeric(as.character(nbeads)))
      )
    }
    
  }) %>% rbindlist()
  
  mean(same_cluster_df$same_cluster)
}

#----------------
#----------------

# Sai's data
cdf <- read.table("../data/sai/clusters.csv", sep = ",", header = TRUE)
bdf <- read.table("../data/sai/possorted_bam.barcodeTranslate.tsv"); colnames(bdf) <- c("Barcode", "bap_id")
bdf$n_beads <- str_sub(bdf$bap_id,-2,-1) %>% as.character() %>% as.numeric
sai_df <- merge(bdf, cdf)

# 10X 5k data
cdf <- read.table("../data/pub_5k/clusters.csv", sep = ",", header = TRUE)
bdf <- read.table("../data/pub_5k/atac_v1_pbmc_5k_possorted_bam.barcodeTranslate.tsv"); colnames(bdf) <- c("Barcode", "bap_id")
bdf$n_beads <- str_sub(bdf$bap_id,-2,-1) %>% as.character() %>% as.numeric
tenx_5k <- merge(bdf, cdf)

# Make a bar plot of similar cluster assignments
if(FALSE){
  prop_cluster_shared_multiplet(sai_df)
  prop_cluster_shared_multiplet(tenx_5k)
  
  permuted_1k_cluster_sai <- sapply(1:100, function(seed){
    prop_cluster_shared_multiplet(sai_df, permute.seed = seed)
  })
  
  permuted_1k_cluster_5k <- sapply(1:100, function(seed){
    prop_cluster_shared_multiplet(tenx_5k, permute.seed = seed)
  })
  
  se <- function(vec){
    return(sqrt(var(vec)))
  }
  
  plot_df <- data.frame(
    what = c("observed", "permuted", "observed", "permuted"),
    data = c("In-house", "In-house", "Public", "Public"),
    mean = c(prop_cluster_shared_multiplet(sai_df), mean(permuted_1k_cluster_sai), 
             prop_cluster_shared_multiplet(tenx_5k), mean(permuted_1k_cluster_5k))*100,
    se = c(0, se(permuted_1k_cluster_sai), 0, se(permuted_1k_cluster_5k)) * 100
  )
  
  p1 <- ggplot(plot_df, aes(x = data, y = mean, fill = what)) +
    geom_bar( color = "black", stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("firebrick", "grey")) +
    pretty_plot(fontsize = 8) + L_border() +
    labs(x = "", y = "% pairs in same cluster") +
    scale_y_continuous(expand = c(0,0),  limits = c(0, 100)) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position = position_dodge(width = 0.7)) 
  cowplot::ggsave(p1, file = "../plots/chromatin_cluster_permuted.pdf", width = 2.5, height = 1.7)
  
}

if(FALSE){
  100-prop_droplets_multiplet(pub_df) * 100
  100-prop_droplets_multiplet(sai_df) * 100
  100-prop_bead_barcodes_in_multiplet(pub_df) * 100
  100-prop_bead_barcodes_in_multiplet(sai_df) * 100
  
  public_df_multiplet_count <- data.frame(prop_bead_barcodes_in_multiplet(pub_df) * 100)
  sai_df_multiplet_count <- data.frame(prop_bead_barcodes_in_multiplet(sai_df) * 100)
  
  p1 <- ggplot(public_df_multiplet_count, aes(x = Var1, y = Freq)) +
    geom_bar(fill = "lightgrey", color = "black", stat = "identity", width = 0.7) +
    pretty_plot(fontsize = 8) + L_border() + ggtitle("Public ") +
    labs(x = "# beads / multiplet", y = "% of total barcodes") +
    scale_y_continuous(expand = c(0,0),  limits = c(0, 87)) 
  
  cowplot::ggsave(p1, file = "../plots/public_barcode_multiplets.pdf", width = 1.6, height = 1.6)
  
  p1 <- ggplot(sai_df_multiplet_count, aes(x = Var1, y = Freq)) +
    geom_bar(fill = "lightgrey", color = "black", stat = "identity", width = 0.7) +
    pretty_plot(fontsize = 8) + L_border() + ggtitle("In-house ") +
    labs(x = "# beads / multiplet", y = "% of total barcodes") +
    scale_y_continuous(expand = c(0,0),  limits = c(0, 87)) 
  
  cowplot::ggsave(p1, file = "../plots/sai_barcode_multiplets.pdf", width = 1.6, height = 1.6)
  
}

if(FALSE){
  rdf <- rbind(hamming_distance_count(sai_df, what = "In-house"),
               hamming_distance_count(tenx_5k, what = "Public"))
  
  ggplot(rdf, aes(x = n_barcodes, y = mean, color = what, group = what)) +
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2) 
  
  pub_df %>% filter(n_beads == 9) %>% arrange(bap_id) %>%
    write.table(row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

