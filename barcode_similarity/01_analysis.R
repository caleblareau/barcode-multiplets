library(data.table)
library(dplyr)
library(BuenColors)

# Function to compute the restricted longest common subsequence
rLCS <- function(a,b){
  
  aa <- strsplit(a, "")[[1]]
  bb <- strsplit(b, "")[[1]]
  
  rle_out <- rle(aa == bb)
  value <- max(rle_out$lengths[rle_out$values])
  value <- ifelse(value >0, value, 0)
  return(value)
}

bcs <- fread("737K-august-2016.txt", header = FALSE)[[1]]

set.seed(1)
n = 1000000
bc_df <- data.frame(
  barcode1 = sample(bcs, size = n, replace = TRUE), 
  barcode2 = sample(bcs, size = n, replace = TRUE), 
  stringsAsFactors = FALSE
)
rLCS_vec <- sapply(1:n, function(i){
  rLCS(bc_df[i,1], bc_df[i,2])
})

pldf <- data.frame(rLCS_vec)
pldf2 <- pldf %>% group_by(rLCS_vec) %>% summarize(density = n()/n)
pldf2$cumSum <- cumsum(pldf2$density)

p1 <- ggplot(pldf2, aes(x = rLCS_vec, y = density*100)) +
  geom_bar(color = "black", fill = "lightgrey", stat = "identity") +
  pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "restricted LCS", y = "% of random pairs") +
  scale_y_continuous() +
  geom_vline(xintercept = 5.5, color = "firebrick")
cowplot::ggsave(p1, file = "null_histogram.pdf", width = 2, height = 2)
