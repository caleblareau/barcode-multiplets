library(dplyr)
library(BuenColors)


p0 <- dpois(0, 1) * 100
p1 <- dpois(1, 1) * 100
p2 <- dpois(2, 1) * 100
p3 <- dpois(3, 1) * 100
p4 = 100 - p0 - p1 - p2 - p3

pois_df <- data.frame(
  n_beads = c("0", "1", "2", "3", "4+"), 
  percentages = c(p0, p1, p2, p3, p4)
)

p1 <- ggplot(pois_df, aes(x = n_beads, y = percentages)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 8) + L_border() +
  ggtitle("Optimal Poisson loading") + labs(x = "# beads / droplet", y = "% of total") +
  scale_y_continuous(expand = c(0,0),  limits = c(0, 66)) 

long_raw <- read.table("FOV10_raw_counts.txt", header = TRUE) %>% reshape2::melt(id.vars = c("numberOfBeads"))

sem_df <- long_raw %>% group_by(variable) %>% 
  mutate(totalBeadsFOV = sum(value)) %>%
  ungroup() %>%
  mutate(prop = value/totalBeadsFOV *100) %>% 
  group_by(numberOfBeads) %>% 
  summarize(valueSum = sum(value), sem = sd(prop)/sqrt(n())) %>% 
  mutate(prop = valueSum/sum(valueSum) * 100)


p2 <- ggplot(sem_df, aes(x = numberOfBeads, y = prop)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 8) + L_border() +
  geom_errorbar(aes(ymin=prop-sem, ymax=prop+sem), width=.2) +
  ggtitle("Empirical loading") + labs(x = "# beads / droplet", y = "% of total") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 66)) 

cowplot::ggsave(cowplot::plot_grid(p1, p2), 
                file = "out_pdfs/barcharts_outPDFs.pdf", width = 3.7, height = 2)

long_raw %>% 
  mutate(b = rep(c(0, 1, 2, 3,4), 10)) %>%
  filter(b != 0) %>% mutate(isOne = b == 1) %>%
  mutate(value2 = b*value) %>%
  group_by(isOne, variable) %>% 
  summarize(totalBeadsFOV = sum(value2)) %>% 
  ungroup() %>%
  group_by(variable) %>% 
  mutate(prop = totalBeadsFOV/sum(totalBeadsFOV) *100) %>% 
  group_by(isOne) %>% 
  summarize(valueSum = sum(totalBeadsFOV), sem = sd(prop)/sqrt(n())) %>% 
  mutate(prop = valueSum/sum(valueSum) * 100) %>% mutate(Affected = c("Yes", "No")) -> odf


p3 <- ggplot(odf, aes(x = Affected, y = prop)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity", width = 0.7) +
  pretty_plot(fontsize = 8) + L_border() +
  geom_errorbar(aes(ymin=prop-sem, ymax=prop+sem), width=.2) +
  ggtitle("  ") + labs(x = "# beads / droplet", y = "% of barcodes") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 75)) 

cowplot::ggsave(p3, 
                file = "out_pdfs/affected_barcodes.pdf", width = 1.5, height = 2)

