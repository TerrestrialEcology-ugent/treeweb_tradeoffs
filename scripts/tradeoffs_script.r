################################################
# R script to reproduce the trade-off analysis #
# 
# Author: Lionel R. Hertzog
################################################

# set wd
setwd("PostDoc_Ghent/Spatialsynthesis_stuff/")

# set a seed
set.seed(20210101)

# load libraries
library(plyr)
library(tidyverse)
library(brms)
library(circlize) # for colorRamp2
library(gridExtra) # for grid.arrange
library(reshape2) # for melt and dcast
library(fmsb) # for radarchart
library(viridis) # for nice colors

# load data
brm_dat <- read.csv("data/tradeoffs_data.csv")

# load helper functions
source("scripts/tradeoffs_function.r")

# fit the model
m_mult <- brm(mvbind(C.stock,Decomp,Biomass,P.germ,Herbivory,Predation,Bird.smi,frass,Sha.veg,Sha.herb,Sha.cara,Sha.spider,Sha.iso,Sha.dip,Sha.bird,Sha.bat)
              ~ speccomb + prox.std + edge.std,
              prior = set_prior("normal(0, 2)", class = "b"), data = brm_dat)


# 1. percent of explained deviance

# grab the names of the variables
nn <- names(brm_dat)[2:17]
nn <- gsub("\\.","",nn)
# reorder the variables
nn <- nn[c(1:4,8,5,7,6,9:16)]
# grab the posterior samples
post <- posterior_samples(m_mult)
# compute the residuals and their sd
residd <- residuals(m_mult, summary = FALSE)
resid_all <- apply(residd, c(1,3), sd)
# apply the get_sd function to all variables
sd_all <- ldply(nn, function(n) get_sd(post, n))

# get the across variables average explained variance
sd_all %>%
  group_by(fraction) %>%
  summarise(estimate = mean(estimate)) -> sd_avg
# nicer variable labels
lbls <- c("C. stock", "Decomposition", "Tree biomass", "Regeneration", "Insect biomass", "Herbivory", "Bird biomass", "Predation",
          "Vegetation div.", "Herbivore div.", "Carabid div.", "Spider div.", "Isopod div.", "Diplopod div.", "Bird div.", "Bat div.")
# the plot
gg_sd <- ggplot(sd_all,aes(x= variable, y = estimate, ymin = conf.low, ymax = conf.high, color = fraction,
                           group = fraction)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_x_discrete( labels = lbls) +
  geom_hline(data = sd_avg, aes(yintercept = estimate, color = fraction), linetype = "dashed") +
  labs(x="", y = "Percentage of explained variance") +
  scale_color_discrete(name = "Effects") +
  theme_classic() 

ggsave("figures/01_explained.pdf", gg_sd, width = 80, height = 70, units = "mm",
       dpi = 600, scale = 2.25)

# 2. trade-offs, Fig. 2

## grab the residual correlation from the fitted model
pp_mult <- posterior_samples(m_mult)

rescor <- pp_mult[grep("rescor",names(pp_mult))] # only keep residual correlation parameters
## data wrangling to get the value in a long data format
rescor %>%
  gather(type, value) %>%
  group_by(type) %>%
  summarise(R = median(value), P = sum(value > 0) / n()) %>% # get median plus "p-value" from the 4000 posterior draws
  mutate(P = ifelse(P < 0.5, 1 - P, P)) %>%
  mutate(Intrisic = ifelse(P > 0.9, "Y","N")) -> m_dd

## fit a null model to get "observed" correlations
m_null <- brm(formula = mvbind(C.stock, Decomp, Biomass, P.germ, Herbivory, Predation, Bird.smi, frass, Sha.veg, Sha.herb, Sha.cara, Sha.spider, Sha.iso, Sha.dip, Sha.bird, Sha.bat) ~
                1, data = brm_dat)

# get the posterior draws
pp_null <- posterior_samples(m_null)
rescor_null <- pp_null[grep("rescor",names(pp_null))]

# some nicer names for the plot
ord <- c("Cstock","Decomp", "Biomass","Pgerm","Herbivory","Predation","Birdsmi","frass",
         "Shaveg","Shaherb","Shacara","Shaspider","Shaiso","Shadip","Shabird","Shabat")

fun <- c("Cstock","Decomp", "Biomass","Pgerm","Herbivory","Predation","Birdsmi","frass")
div <- c("Shaveg","Shaherb","Shacara","Shaspider","Shaiso","Shadip","Shabird","Shabat")


## data wrangling to put residual correlation estimate in plotting format
rescor_null %>%
  gather(type, value) %>%
  group_by(type) %>%
  summarise(R = median(value), P = sum(value > 0) / n()) %>% # compute observed median correlation and "p-value"
  mutate(P = ifelse(P < 0.5, 1 - P, P)) %>%
  mutate(Present = ifelse(P > 0.9, "Y","N")) %>% # if p-value higher than 0.9 flag the correlation
  left_join(m_dd[,c(1,4)], by = "type") %>% # add residual correlation after fitting the model
  mutate(Res = paste(Present, Intrisic, sep = "_")) %>% # put together the before and after model fitting residual correlation flags
  separate(type, c(NA,"From","To"), sep  = "__") %>%
  mutate(From = factor(From, levels = ord), To = factor(To, levels = ord)) %>% # re-ordering of the labels
  mutate(label = ifelse(Res == "Y_Y",sprintf("underline(%s)",round(R,2)),round(R,2))) %>% # if intrisic correlation underline
  mutate(label = ifelse(Res == "Y_N", sprintf("italic('%s')",round(R,2)), label)) %>% # if extrinsic correlation italized
  mutate(size = ifelse(Res %in% c("Y_N","Y_Y"), 2, 1)) %>% # bigger plotting size if correlation present
  mutate(type_f = ifelse(From %in% fun, "function","diversity"),
         type_t = ifelse(To %in% fun, "function","diversity")) %>%
  unite("type",c(type_f,type_t),sep="_") %>%
  mutate(color = colorRamp2(c(-0.55,0,0.55),c("red","white","blue"))(R)) %>%
  mutate(From = revalue(From, c(Cstock = "C stock", Decomp = "Decomposition", 
                                Biomass = "Tree biomass", Pgerm = "Regeneration", Birdsmi = "Bird biomass",
                                frass = "Insect biomass",
                                Shaveg = "Vegetation", Shaherb = "Herbivore", Shacara = "Carabid", Shaspider = "Spider", 
                                Shaiso = "Isopod", Shadip = "Diplopod", Shabird = "Bird", Shabat = "Bat"))) %>% # re-name variables for better plotting visibility
  mutate(To = revalue(To, c(Cstock = "C stock", Decomp = "Decomposition", Biomass = "Tree biomass", Pgerm = "Regeneration", Birdsmi = "Bird biomass", frass = "Insect biomass",
                            Shaveg = "Vegetation", Shaherb = "Herbivore", Shacara = "Carabid", Shaspider = "Spider", Shaiso = "Isopod", Shadip = "Diplopod", Shabird = "Bird", Shabat = "Bat"))) %>%
  mutate(type = revalue(type, c(diversity_diversity = "a) Diversity - Diversity", function_function = "b) Function - Function", function_diversity = "c) Function - Diversity"))) -> obs_dd

## the plot, note that the values may slightly differ from the published version due to stochastic sampling in the model
gg_cor <- ggplot(obs_dd,aes(x=From,y=To,fill=color)) +
  geom_tile() +
  geom_text(aes(label = label, size = size), parse = TRUE, show.legend = FALSE) +
  scale_fill_identity() +
  facet_wrap(~type,scales = "free", ncol = 3) +
  scale_size(range = c(2,4)) +
  labs(x = "", y = "") +
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), axis.ticks = element_blank()) 


ggsave("figures/02_tradeoffs.pdf",gg_cor,width=180, height = 60, units = "mm", dpi = 600,
       scale = 1.75)


# 3. Desirability, Fig. 3

## read-in importance score, note you can create your own !!
imp <- read.csv("data/importance_score.csv")

## 3.1 desirability at plot-scale
newdat_plot <- expand.grid(speccomb = c("fsyl","qrob","qrub","all"), prox.std = c(0, -1.37), edge.std = c(2, -1.1))[c(1:4,13:16),] # data to predict on

## get model predictions
pred_plot <- posterior_linpred(m_mult, newdata = newdat_plot)
## get desirability scores from productivist and conservationist perspective
dd_plotf <- desirab(pred_plot, imp$direction_f, imp$weight_f)
dd_plotf$perspective <- "productivist"
dd_plotc <- desirab(pred_plot, imp$direction_c, imp$weight_c)
dd_plotc$perspective <- "conservationist"
## put the two perspective together
dd_plot <- rbind(dd_plotf, dd_plotc)
dd_plot$fragmentation <- rep(rep(c("high frag.","low frag."),each = 4), 4) # some labelling for the plot
dd_plot$speccomb <- rep(newdat_plot$speccomb, 4)
dd_plot$X <- rep(1:4, 8) + c(rep(rep(c(-0.15, 0.05), each = 4), 2),rep(rep(c(-0.05, 0.15), each = 4), 2)) # x coords on the plot
dd_plot$top <- "plot-level"
# the plot
gg_plot <- ggplot(dd_plot, aes(x = X, y = Median, ymin = LCI, ymax = UCI, color = fragmentation, shape = perspective)) +
  geom_linerange() +
  geom_point(size = 2.5) +
  facet_grid(type ~ top) +
  scale_x_continuous(breaks = 1:4, labels = c("fsyl","qrob","qrub","all")) +
  labs(x = "", y = "Desirability scores (with 80% CrI)") +
  theme(legend.position = "none", strip.text.y = element_blank()) 

## 3.2 Desirability at landscape-scale

## generate hypothetical landscapes
newdat_land1 <- data.frame(speccomb = rep(c("fsyl","qrob","qrub"), times = c(17,18,18)), edge.std = 2, prox.std = 0) # monoculture / high frag.
newdat_land2 <- data.frame(speccomb = rep("all", 53), edge.std = 2, prox.std = 0, dens.std = 0) # mixture / high frag
newdat_land3 <- data.frame(speccomb = rep(c("fsyl","qrob","qrub"), times = c(17,18,18)), edge.std = -1.11, prox.std = -1.37, dens.std = 0) # monoculture / low frag
newdat_land4 <- data.frame(speccomb = rep("all", 53), edge.std = -1.11, prox.std = -1.37, dens.std = 0) # mixture / low frag

## get the model predictions
pred_1 <- posterior_linpred(m_mult, newdata = newdat_land1)
pred_2 <- posterior_linpred(m_mult, newdata = newdat_land2)
pred_3 <- posterior_linpred(m_mult, newdata = newdat_land3)
pred_4 <- posterior_linpred(m_mult, newdata = newdat_land4)
pred_l <- list("monoculture\nhigh frag."=pred_1, "mixtures\nhigh frag." = pred_2, "monoculture\nlow frag." = pred_3, "mixtures\nlow frag." = pred_4)

## compute desirability
dd_f <- ldply(pred_l, function(x) desirab(x, imp$direction_f, imp$weight_f, sum = TRUE)) # for productivist
dd_c <- ldply(pred_l, function(x) desirab(x, imp$direction_c, imp$weight_c, sum = TRUE)) # for conservationist
dd_f$perspective <- "productivist"
dd_c$perspective <- "conservationist"
# combine the two perspectives
dd_all <- rbind(dd_f,dd_c)
## some cleaning up before plot
dd_all <- separate(dd_all, ".id", c("tree_composition", "fragmentation"), sep = "\n")
dd_all$X <- rep(1:2, each = 2) + rep(c(-0.15, 0.05, -0.05,0.15), each = 4)
dd_all$top <- "landscape"

gg_land <- ggplot(dd_all,aes(x=X,y=Median, ymin = LCI, ymax = UCI, color = fragmentation, shape = perspective)) +
  geom_linerange() +
  geom_point(size = 2.5) +
  facet_grid(type ~ top) +
  labs(x = "", y = "") +
  scale_x_continuous(breaks = 1:2, labels = c("monoculture","mixture"), limits = c(0.5,2.5))

gg_desirab <- grid.arrange(gg_plot, gg_land, bottom = "Tree species composition", widths = c(6.5, 8.5))

ggsave("figures/03_desirab2.png", gg_desirab, width = 7.7, height = 6)

