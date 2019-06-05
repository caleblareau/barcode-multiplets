library(dplyr)
library(BuenColors)
library(compoisson)

getThree_fromZero <- function(p_zero){
  lambda <- -log(p_zero)
  p_one <- lambda * exp(-1*lambda)
  p_more <- 1 - p_zero - p_one
  return(c(p_zero, p_one, p_more))
}

sapply(seq(0.001, 0.999, 0.0001), getThree_fromZero) %>% t() %>% data.frame() -> odf
colnames(odf) <- c("p0", "p1", "pM")

p1 <- ggplot(odf, aes(x = pM*100,y = p1*100, color = p0*100)) + 
  geom_point() + pretty_plot() + L_border() +
  labs(x = "% drops with multiple beads", y = "% drops with 1 bead", color = "% with 0 beads") +
  theme(legend.position = "bottom")
cowplot::ggsave(p1, file = "poisson_beads.pdf", width = 3, height = 3)

# Get MLE from Poisson
cmat <- matrix(data = c(0, 13.5, 1, 81.5,2,5), ncol = 2, byrow = TRUE)
fit <- com.fit(cmat)
p_zero <- dcom(0, fit$lambda, fit$nu)
com.mean(fit$lambda, fit$nu)
com.var(fit$lambda, fit$nu)


# Seed the variance from the MLE above
COM_tenX <- function(mean_desired){
  variance_fixed <- 0.1779
  find_two <- function(l_n){
    lambda <- l_n[1]
    nu <- l_n[2]
    if(nu < 0 | lambda < 0) return(9999)
    abs(com.var(lambda, nu) - variance_fixed) + abs(com.mean(lambda, nu) - mean_desired)
  }
  
  values <- optim(c(6,6), find_two)$par
  p_zero <- dcom(0, lambda = values[1], nu = values[2])
  p_one <- dcom(1, lambda = values[1], nu = values[2])
  p_more <- 1-p_zero - p_one
  return(c(p_zero, p_one, p_more))
}
sapply(c(seq(0.1, 1.5, 0.02), seq(1.5, 3, 0.1)), COM_tenX) %>% t() %>% data.frame() -> odf

colnames(odf) <- c("p0", "p1", "pM")
p1 <- ggplot(odf %>% filter(p1 > 0.2 | pM > 0.25), aes(x = pM*100,y = p1*100, color = p0*100)) + 
  geom_point() + pretty_plot() + L_border() +
  labs(x = "% drops with multiple beads", y = "% drops with 1 bead", color = "% with 0 beads") +
  theme(legend.position = "bottom")
cowplot::ggsave(p1, file = "10X_beads.pdf", width = 3, height = 3)
