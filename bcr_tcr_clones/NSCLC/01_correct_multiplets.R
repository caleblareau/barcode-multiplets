library(SummarizedExperiment)
library(dplyr)
library(reshape2)
library(data.table)
library(BuenColors)
"%ni%" <- Negate("%in%")

one_multinom_generator <- function(prob_vec){
  which(1==rmultinom(n = 1, size = 1,prob = prob_vec)[,1]) %>% unname()
}

make_plot_corrected_clonotypes <- function(file, what){
  
  # Filter down for consensus clonotypes
  clones <- read.table(file,
                       sep = ",", header = TRUE, stringsAsFactors = FALSE)
  clones %>% filter(is_cell == "True"  & raw_clonotype_id != "None") %>%
    dplyr::select(barcode, raw_clonotype_id) %>% distinct() %>%
    group_by(raw_clonotype_id) %>%
    mutate(n = n()) %>% arrange(raw_clonotype_id) -> clone_df
  
  clone_df%>% group_by(n) %>% summarize(clone_count = n()) %>% mutate(uniq_clone_count = clone_count/n) -> df
  clone_df %>% dplyr::select(c("raw_clonotype_id", "n")) %>% unique() %>% data.frame() -> unique_clone_df
  
  total_barcodes <- sum(df$clone_count)
  
  
  classify <- function(df){
    df$group <- cut(df$n, breaks = c(1, 2, 5, 10, 100), right = FALSE)
    df %>% group_by(group) %>%
      summarize(count = sum(clone_count)) %>%
      mutate(prop = count / sum(count) * 100)
  }
  
  simulate_bead_doublet_collapse <- function(seed_exp){
    
    # Generate putative multiple draw
    set.seed(seed_exp)
    probs <- c(0.85, 0.1, 0.02, 0.02, 0.01)
    multiplet_draw <- rmultinom(n = total_barcodes, size = 1,prob = probs) %>% melt() %>%
      filter(value != 0) %>% pull(Var1) 
    
    # Remove 1s that we know are real and recompute probabilities
    ones_idx <- head(which(multiplet_draw==1), df %>% filter(n == 1) %>% pull(clone_count))
    multiplet_draw_singlets_removed <- sample(multiplet_draw[1:length(multiplet_draw) %ni% ones_idx])
    probs_adjusted <- table(sort(multiplet_draw_singlets_removed))/length(multiplet_draw_singlets_removed) %>% as.numeric()
    
    # Loop through each clone and compute number of droplets / cells needed to make it up
    unique_clone_df$corrected_count <- sapply(1:dim(unique_clone_df)[1], function(j){
      total <- unique_clone_df[j,2] %>% as.numeric()
      if(total == 1){
        1
      } else {
        new_total <- 0
        drawn = 0
        while(new_total < total){
          obs <- one_multinom_generator(probs_adjusted)
          new_total = new_total + obs
          drawn = drawn + 1
        }
        drawn
      }
      
    })
    
    unique_clone_df%>% group_by(corrected_count) %>% summarize(clone_count = n()) %>% 
      mutate(n = corrected_count) %>%
      mutate(clone_count = corrected_count*clone_count) %>%
      mutate(uniq_clone_count = clone_count/n) %>% classify() -> corrected_df
    
    clone_FDR <- sum(unique_clone_df$n > 1 & unique_clone_df$corrected_count == 1)/sum(unique_clone_df$n > 1)
    
    corrected_df$seed <- seed_exp
    l2 <- list(corrected_df, clone_FDR)
    return(l2)
  }
  
  # Compare against original
  lists_o_lists <- lapply(1:100, simulate_bead_doublet_collapse) 
  
  # Aggregate over simulation
  lapply(lists_o_lists, function(x){
    x[[1]]
  }) %>%
    rbindlist() %>% data.frame() %>% group_by(group) %>%
    summarize(mean = mean(prop), se = sqrt(var(prop)/10)) %>%
    mutate(prop = mean, count = 1, what = "adjusted") -> adjusted_df
  
  #Aggregate for CLone FDR
  
  clone_fdr <- sapply(lists_o_lists, function(x){
    x[[2]]
  }) %>% mean()
  print(paste0(what, " clone FDR: ", as.character(clone_fdr)))
  
  # Reaggregate
  clone_df%>% group_by(n) %>% summarize(clone_count = n()) %>% mutate(uniq_clone_count = clone_count/n) %>% classify()  %>%
    mutate(se = 0, what = "observed", mean = 0) -> original_df
  
  total_df <- rbind(adjusted_df, original_df)
  total_df$what <- factor(total_df$what, levels = c("observed", "adjusted"))
  
  p1 <- ggplot(total_df, aes(x = group, y = prop, fill = what)) +
    geom_bar(color = "black", stat = "identity", width = 0.5,position = "dodge") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position = position_dodge(width = 0.5)) +
    pretty_plot(fontsize = 7) + L_border() +
    scale_y_continuous(expand = c(0,0), limits = c(70, 90)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "# Cells / Clone", y = "% of Cells in Library") +
    scale_fill_manual(values = c("darkgrey", "lightgrey"))
  
  
  p2 <- ggplot(total_df, aes(x = group, y = prop, fill = what)) +
    geom_bar(color = "black", stat = "identity", width = 0.5,position = "dodge") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position = position_dodge(width = 0.5)) +
    pretty_plot(fontsize = 7) + L_border() +
    scale_y_continuous(expand = c(0,0), limits = c(0,20)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "# Cells / Clone", y = "% of Cells in Library") +
    scale_fill_manual(values = c("darkgrey", "lightgrey")) +
    theme(legend.position = "none")
  
  
  cowplot::ggsave(cowplot::plot_grid(p2, ncol = 1), file = paste0("plots/",what,"_adjusted.pdf"), width = 2, height = 1.3)
  what
}

make_plot_corrected_clonotypes("../data/nsclc/vdj_v1_hs_nsclc_b_all_contig_annotations.csv", "BCR")
make_plot_corrected_clonotypes("../data/nsclc/vdj_v1_hs_nsclc_t_all_contig_annotations.csv", "TCR")
