library(dplyr)
library(data.table)
library(stringr)
library(BuenColors)

import_rate <- function(i){
  file <- paste0("../data/pub_5k/diff_topk/pbmc5k-top",i,"k.barcodeTranslate.tsv")
  ns <- str_split_fixed(fread(file, header = FALSE)[[2]], "_", 3)[,3] 
  sum(ns != "N01")/length(ns)
}

plot_df <- data.frame(
  multiplet_rate = c(sapply(5:10, import_rate),0.178 )*100,
  n_barcodes = c(5:10, 5.2)
)

p1 <- ggplot(plot_df, aes(x = n_barcodes, y = multiplet_rate)) +
  geom_point() + geom_line() +
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(limits = c(0, 30)) +
  labs(x = "Number of barcodes (x 103)", y = "Barcode multiplets (%)")
cowplot::ggsave(p1, file = "../plots/multiplet_rate_plot.pdf", width = 2, height =2 )
