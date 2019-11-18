library(dplyr)
library(BuenColors)

df <- data.frame(
  technology = c("v1", "v11", "v1", "v11"), 
  what = c("zcomplex", "zcomplex", "multi", "multi"), 
  rate = c(83.5, 86.1, 16.5, 13.9)
)

p1 <- ggplot(df, aes(x = technology, y = rate, fill = what)) +
  geom_bar(color = "black", stat = "identity") +
  pretty_plot(fontsize = 8) + L_border() +
  scale_fill_manual(values = c("darkgrey", "lightgrey"))+
  labs(x = "Chip", y = "% of total multiplets") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100))  +
  theme(legend.position = "none")

cowplot::ggsave(p1, file = "plot_stacked.pdf", width = 1.5, height = 2)
