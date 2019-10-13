library(dplyr)
library(BuenColors)
library(ggbeeswarm)

long_raw1 <- read.table("FOV10_raw_counts.txt", header = TRUE) %>% reshape2::melt(id.vars = c("numberOfBeads"))
long_raw2 <- read.table("FOV10_raw_counts2.txt", header = TRUE) %>% reshape2::melt(id.vars = c("numberOfBeads"))
long_raw3 <- read.table("FOV10_raw_counts3.txt", header = TRUE) %>% reshape2::melt(id.vars = c("numberOfBeads"))

make_sem_df <- function(long_raw){
  sem_df <- long_raw %>% group_by(variable) %>% 
    mutate(totalBeadsFOV = sum(value)) %>%
    ungroup() %>%
    mutate(prop = value/totalBeadsFOV *100) %>% 
    group_by(numberOfBeads) %>% 
    summarize(valueSum = sum(value), sem = sd(prop)/sqrt(n())) %>% 
    mutate(prop = valueSum/sum(valueSum) * 100)
  sem_df
}

make_odf <- function(long_raw){
  
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
  odf
}

sem_df1 <- make_sem_df(long_raw1)
sem_df2 <- make_sem_df(long_raw2)
sem_df3 <- make_sem_df(long_raw3)

odf1 <- make_odf(long_raw1)
odf2 <- make_odf(long_raw2)
odf3 <- make_odf(long_raw3)

p1 <- ggplot(sem_df1, aes(x = numberOfBeads, y = prop)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 8) + L_border() +
  geom_errorbar(aes(ymin=prop-sem, ymax=prop+sem), width=.2) +
  ggtitle("Rep 1") + labs(x = "# beads / droplet", y = "% of total") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 95)) 

p2 <- ggplot(sem_df2, aes(x = numberOfBeads, y = prop)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 8) + L_border() +
  geom_errorbar(aes(ymin=prop-sem, ymax=prop+sem), width=.2) +
  ggtitle("Rep 2") + labs(x = "# beads / droplet", y = "% of total") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 95)) 

p3 <- ggplot(sem_df3, aes(x = numberOfBeads, y = prop)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 8) + L_border() +
  geom_errorbar(aes(ymin=prop-sem, ymax=prop+sem), width=.2) +
  ggtitle("Rep 3 (Regev)") + labs(x = "# beads / droplet", y = "% of total") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 95)) 

cowplot::ggsave(cowplot::plot_grid(p1, p2, p3, nrow = 1), 
                file = "out_pdfs/barcharts_outPDFs.pdf", width = 5.5, height = 2)

plot_df <- data.frame(
  dataset = rep(c("Rep1", "Rep2", "Rep3", "Zheng"), 2),
  value = c(27.1, 4.5, 2.5, 11.7, 72.9, 95.5, 97.5, 88.3),
  doublet = c(rep("Yes", 4), rep("No", 4))
)
mean_df <- plot_df %>% group_by(doublet) %>% summarize(sem = sd(value)/sqrt(n()), value = mean(value))
pout <- ggplot(plot_df, aes(x = doublet, y = value, color = dataset)) +
  geom_bar(data = mean_df, color = "black", fill = "lightgrey", stat = "identity", width = 0.7)+
  geom_errorbar(data = mean_df, aes(ymin=value-sem, ymax=value+sem), width=.2, color = "black") +
  geom_quasirandom(size = 0.5) + 
  pretty_plot(fontsize = 8) + L_border() +
  ggtitle("  ") + labs(x = "affected by doublet", y = "% of barcodes", color = "") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) 

cowplot::ggsave(pout, 
                file = "out_pdfs/affected_barcodes.pdf", width = 2.3, height = 2)

beads_df <- data.frame(
  dataset = rep(c("Rep1", "Rep2", "Rep3", "Zheng"), 3),
  value = c(28.1, 14.3, 7.11, 15, 63.5, 84.2, 92.1, 80, 8.36, 1.5, 0.8, 5),
  nBeads = c(rep("0", 4), rep("1", 4), rep("2+", 4))
)

mean_df2 <- beads_df %>% group_by(nBeads) %>% summarize(sem = sd(value)/sqrt(n()), value = mean(value))

pouts2 <- ggplot(beads_df, aes(x = nBeads, y = value, color = dataset)) +
  geom_bar(data = mean_df2, color = "black", fill = "lightgrey", stat = "identity", width = 0.7)+
  geom_errorbar(data = mean_df2, aes(ymin=value-sem, ymax=value+sem), width=.2, color = "black") +
  geom_quasirandom(size = 0.5) + 
  pretty_plot(fontsize = 8) + L_border() +
  ggtitle("  ") + labs(x = "# beads / droplet", y = "% of droplets", color = "") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) + theme(legend.position = "none")

cowplot::ggsave(pouts2, 
                file = "out_pdfs/beads_per_droplet.pdf", width = 2, height = 2)


