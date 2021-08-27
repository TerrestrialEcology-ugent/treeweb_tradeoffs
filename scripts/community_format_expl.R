# script to explore and format community data to use in synthesis II

## set wd
setwd("~/PostDoc_Ghent/Spatialsynthesis_stuff/")

## load libraries
library(tidyverse)

## load data
bird <- read.csv("data/treeweb_bird.csv")
bat <- read.csv("data/treeweb_bat.csv")
spider <- read.csv("data/treeweb_spider_ground.csv")
miner <- read.csv("data/treeweb_miner.csv")
gall <- read.csv("data/treeweb_gall.csv")
veg <- read.csv("data/treeweb_vegetation.csv")
cara <- read.csv("data/treeweb_carabids.csv")
iso <- read.csv("data/treeweb_isopods_all.csv",stringsAsFactors = FALSE)
iso_tmp <- read.csv("data/treeweb_isopodes_may_tmp.csv",stringsAsFactors = FALSE)
dip <- read.csv("data/treeweb_diplopods.csv")

## a function for some div measure derivation
abg <- function(comm, name = "A"){ # comm needs to be in relative values in cells
  N <- nrow(comm)
  alpha <- exp(sum(apply(comm,1,function(x) -(1/N) * sum(x * log(x),na.rm=TRUE))))
  gamma <- exp(sum(apply(comm,2,function(x) -sum(x * (1/N)) * log(sum(x * (1/N))))))
  beta <- gamma / alpha
  turnover <- (beta - 1) / (N - 1)
  return(data.frame(group = name, alpha = alpha, beta = beta, gamma = gamma, turnover = turnover))
}

## exploration and formatting of bird data
bird %>%
  rename(id_plot = plot) %>%
  mutate(id_plot = as.character(id_plot)) %>%
  group_by(id_plot) %>%
  mutate(rel_abun = count / sum(count)) %>%
  summarise(Rich = sum(count > 0), 
            Sha = exp(-sum(rel_abun * log(rel_abun), na.rm=TRUE)),
            Simp = sum(rel_abun ** 2) ** (-1),
            Abun = sum(count)) -> bird_dd

cor(bird_dd[,-1]) # high correlation between div measures
hist(bird_dd$Rich)
hist(bird_dd$Sha)
hist(bird_dd$Simp)
hist(bird_dd$Abun)

## exploration and formatting of bat data
bat %>%
  mutate(id_plot = as.character(id_plot)) %>%
  group_by(id_plot) %>%
  mutate(rel_abun = Calls / sum(Calls)) %>%
  summarise(Rich = sum(Calls > 0), 
            Sha = exp(-sum(rel_abun * log(rel_abun), na.rm=TRUE)),
            Simp = sum(rel_abun ** 2) ** (-1),
            Abun = sum(Calls)) -> bat_dd

cor(bat_dd[,-1]) # low correlation between richness and the other div measures
hist(bat_dd$Rich)
hist(bat_dd$Sha) # un-even communities
hist(bat_dd$Simp)
hist(bat_dd$Abun)

## exploration and formatting of spider data
spider %>%
  gather("species","count",-X) %>%
  rename(id_plot = X) %>%
  mutate(id_plot = as.character(id_plot)) %>%
  group_by(id_plot) %>%
  mutate(rel_abun = count / sum(count)) %>%
  summarise(Rich = sum(count > 0), 
           Sha = exp(-sum(rel_abun * log(rel_abun), na.rm=TRUE)),
           Simp = sum(rel_abun ** 2) ** (-1),
           Abun = sum(count)) -> spider_dd

cor(spider_dd[,-1]) # again low correlation richness - shannon
hist(spider_dd$Rich)
hist(spider_dd$Sha)
hist(spider_dd$Simp)
hist(spider_dd$Abun)

## exploration and formatting of leaf herbivore data
miner %>%
  group_by(id_plot,spec) %>%
  summarise(abun = sum(abun)) -> miner_dd
miner_dd$spec <- as.character(miner_dd$spec)

gall %>%
  group_by(id_plot, spec) %>%
  summarise(abun = sum(abun)) -> gall_dd
gall_dd$spec <- as.character(gall_dd$spec)

herb <- rbind(miner_dd, gall_dd)
herb %>%
  ungroup() %>%
  mutate(id_plot = as.character(id_plot)) %>%
  group_by(id_plot) %>%
  mutate(rel_abun = abun / sum(abun)) %>%
  summarise(Rich = sum(abun > 0), 
            Sha = exp(-sum(rel_abun * log(rel_abun), na.rm=TRUE)),
            Simp = sum(rel_abun ** 2) ** (-1),
            Abun = sum(abun)) -> herb_dd

cor(herb_dd[,-1])
hist(herb_dd$Abun)
hist(herb_dd$Rich)
hist(herb_dd$Sha)
hist(herb_dd$Simp)

## exploration and formatting of understorey vegetation data
veg %>%
  mutate(id_plot = as.character(id_plot)) %>%
  filter(layer == "herb") %>%
  group_by(id_plot, species_abbr) %>%
  summarise(cover = sum(cover)) %>%
  group_by(id_plot) %>%
  mutate(rel_abun = cover / sum(cover)) %>%
  summarise(Rich = sum(cover > 0), 
            Sha = exp(-sum(rel_abun * log(rel_abun), na.rm=TRUE)),
            Simp = sum(rel_abun ** 2) ** (-1),
            Abun = sum(cover) / 5) -> veg_dd # divided by 5 to account for the 5 subplots

cor(veg_dd[,-1])
hist(veg_dd$Abun)
hist(veg_dd$Rich)
hist(veg_dd$Sha)
hist(veg_dd$Simp)

# epxloration and formatting of carabid data
cara %>%
  filter(CODE.Pallieter != "eentje zonder naam") %>%
  filter(Species != "geen carabidae") %>%
  separate(CODE.Pallieter, c("id_plot","Month","replicate")) %>%
  group_by(id_plot, Species) %>%
  summarise(Abun = sum(Aantal)) %>%
  group_by(id_plot) %>%
  mutate(rel_abun = Abun / sum(Abun)) %>%
  summarise(Rich = n(),
            Sha = exp(-sum(rel_abun * log(rel_abun),na.rm=TRUE)),
            Simp = sum(rel_abun ** 2, na.rm=TRUE) ** (-1),
            Abun = sum(Abun, na.rm=TRUE)) -> cara_dd

# exploration and formatting of isopodes
# some formatting of names
iso$Genus <- paste0(substring(iso$Genus,1,1),tolower(substring(iso$Genus,2,nchar(iso$Genus))))
iso$Species <- tolower(iso$Species)

# put back with the wrongly placed individuals
iso_tmp %>%
  separate(Sample.ID,c("Plot","Date","Replicate","Drop")) -> iso_tmp
  
isot <- iso_tmp[-9,]
isot <- isot[,-4]
names(isot)[6:7] <- c("Amount","Sex")
names(iso)[5] <- "species"

iso <- rbind(iso,isot)

iso %>%
  rename(id_plot = Plot) %>%
  mutate(id_plot = as.character(id_plot)) %>%
  unite(Genus_species,Genus,species) %>%
  group_by(id_plot,Genus_species) %>%
  summarise(Abun = sum(Amount)) %>%
  group_by(id_plot) %>%
  mutate(rel_abun = Abun / sum(Abun)) %>%
  summarise(Rich = n(),
            Sha = exp(-sum(rel_abun * log(rel_abun),na.rm=TRUE)),
            Simp = sum(rel_abun ** 2, na.rm=TRUE) ** (-1),
            Abun = sum(Abun, na.rm=TRUE)) -> iso_dd

# iso_dd <- rbind(iso_dd,data.frame(id_plot="11",Rich=0,Sha=0,Simp=0,Abun=0))

# exploring and formatting of diplopodes
dip %>%
  separate(Sample.ID,c("id_plot","Month","Replicate","Drop")) %>%
  unite(Genus_species,Genus,species) %>%
  group_by(id_plot, Genus_species) %>%
  summarise(Abun = sum(Count)) %>%
  group_by(id_plot) %>%
  mutate(rel_abun = Abun / sum(Abun)) %>%
  summarise(Rich = n(),
            Sha = exp(-sum(rel_abun * log(rel_abun),na.rm=TRUE)),
            Simp = sum(rel_abun ** 2, na.rm=TRUE) ** (-1),
            Abun = sum(Abun, na.rm=TRUE)) -> dip_dd

dip_dd <- rbind(dip_dd,data.frame(id_plot=c("11","38"),Rich=0,Sha=0,Simp=0,Abun=0))


# put all together
list(veg_dd,herb_dd,cara_dd,spider_dd,iso_dd,dip_dd,bird_dd,bat_dd) %>%
  reduce(left_join, by = "id_plot") -> div_all

names(div_all)[-1] <- paste(rep(c("Rich","Sha","Simp","Abun"),rep=8),rep(c("veg","herb","cara","spider","iso","dip","bird","bat"),each=4),sep="_")



# save that file
write.table(div_all,"data/community_data_formatted.csv",row.names=FALSE,sep=",")


### dev, some cutting into categories
pred_dat %>%
  mutate(prox_cut = cut_number(frag_prox,2), edge_cut = cut_number(plot_edge100m,2)) -> dd

# compute prop of tit
tits <- c("Blue Tit - Cyanistes caeruleus", "Great Tit - Parus major")
bird %>%
  filter(place == "inside") %>%
  group_by(plot) %>%
  summarise(tot = sum(count),
            tit = sum(count[species %in% tits])) %>%
  mutate(prop = tit / tot) -> dd
