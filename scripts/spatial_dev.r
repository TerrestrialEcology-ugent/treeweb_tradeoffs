# some spatial musing

setwd("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/")

# load data
# dat_std <- read.csv("../Data/synthesis_responsedata_std.csv",sep=" ")
# dat <- read.csv("../Data/synthesis_responsedata_raw.csv",sep=" ")
# div_dd <- read.csv("../Data/synthesis_expldata_raw.csv",sep=" ")
# obs <- read.csv("../Data/Plot_diversity_20170113.csv")

# create neighbor matrix, two plots in same fragments are neighbor
W <- array(0,c(53,53))
distance <- as.matrix(dist(div_dd[,c(14,14)]))
W[distance == 0] <- 1

dat <- data.frame(biom = dat_std$Biomass, rich = div_dd$specrich)

# fit a simple model
fit <- brm(biom ~ rich, data = dat, autocor = cor_car(W),iter = 3000,control = list(adapt_delta = 0.9))


exp_quad <- function(xy, a = 1, p = 0.5, sigma = 1){
  dd <- as.matrix(dist(xy))
  K <- a * exp(-(1/(2*p)) * dd)
  diag(K) <- sigma
  plot(as.numeric(dd),as.numeric(K))
}


# some Stan models
n <- 100
xy <- data.frame(x=runif(n,-2, 2),y=runif(n,-2,2))
a_sq <- 1
p_sq <- 1
s_sq <- 1

# m <- stan("~/Documents/PostDoc_Ghent/STAN_stuff/Models/gp_simu.stan", data = list(N=n,x=x,D=D,a_sq=a_sq,p_sq=p_sq,s_sq=s_sq),iter = 100)

# function to generate cov_exp_quad

cov_exp_quad <- function(xy, alpha, rho, sigma){
  K <- alpha ** 2 * exp((-(1 / 2 * rho **2)) * as.matrix(dist(xy)) ** 2)
  diag(K) <- sigma
  return(K)
}

# a function to plot spatial covariance as a function of distance 
cov_dist <- function(dist, alpha, rho){
  cov_est <- alpha ** 2 * exp((-(1 / 2 * rho **2)) * (dist ** 2))
  return(cov_est)
}


K <- cov_exp_quad(xy, a_sq, p_sq, s_sq)

# generate a sample
y <- mvrnorm(mu = rep(0, n), Sigma = K)

# fit the model
m <- stan("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/model/gp_2D_basic.stan", data = list(N = n, D = 2, x = xy, y = y)) # works!

# plot the dist - covariance relation
post_C <- extract(m, c("alpha","rho"))
PP <- sapply(sample(1:4000,100,replace=FALSE), function(iter) cov_dist(seq(1,10,length.out=10),post_C$alpha[iter],post_C$rho[iter]))
# get the prediction for the median
PP_med <- cov_dist(seq(1,10,length.out = 10), median(post_C$alpha), median(post_C$rho))

plot(seq(1,10,length.out = 10), PP_med, lwd = 2, type = "l", ylim = c(0, max(PP)))
for(i in 1:100){
  lines(seq(1,10, length = 10), PP[,i],col = "grey10", lwd = 0.5)
}


# derive posterior predicted correlations
post_K <- extract(m, "K")
# posterior median
K_med <- apply(post_m$K, c(2,3), median)
# turn to corr matrix
R_med <- cov2cor(K_med)

# plot the network
plot(y ~ x, xy, pch = 16, col = "blue")
for(i in 1:100){
  for(j in 1:100){
    if(i < j)
      lines(c(xy$x[i],xy$x[j]),c(xy$y[i],xy$y[j]),col = col.alpha("black",R_med[i, j]))
  }
}

# try with some predictors
X <- matrix(rep(1,n),ncol=1)
beta <- 2

y <- mvrnorm(mu = X * beta, Sigma = K)

m <- stan("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/model/gp_2D_pred.stan",
          data = list(N = n, D = 2, P = ncol(X), coords = xy, X = X, y = y))

# try with some real data
func <- read.csv("data/spatsynth_fun.csv")
expl <- read.csv("data/spatsynth_expl.csv")
coords <- read.csv("data/coords.csv")

# scale edge eff
expl$edge_s <- scale(expl$edge)

# a first model
X <- model.matrix(~ specrich * edge_s, expl)

m_1 <- stan("model/gp_2D_pred_v2.stan",
            data = list(N=53,D=2,P=ncol(X),coords=coords[,-1],X=X,y=func$Predation),
            control = list(max_treedepth = 15))

# plot predictors
plot(m_1, pars = "beta")
plot(m_1, pars = c("rho","alpha","sigma"))

# plot the dist - covariance relation
post_C <- extract(m_1, c("alpha","rho"))
PP <- sapply(sample(1:4000,100,replace=FALSE), function(iter) cov_dist(seq(1,100,length.out=100),post_C$alpha[iter],post_C$rho[iter]))
# get the prediction for the median
PP_med <- cov_dist(seq(1,100,length.out = 100), median(post_C$alpha), median(post_C$rho))

plot(seq(1,100,length.out = 100), PP_med, lwd = 2, type = "l", ylim = c(0, max(PP)))
for(i in 1:100){
  lines(seq(1,100, length = 100), PP[,i],col = "grey10", lwd = 0.5)
}

# derive posterior predicted correlations
post_K <- extract(m_1, "K")
# posterior median
K_med <- apply(post_K$K, c(2,3), median)
# turn to corr matrix
R_med <- cov2cor(K_med)

# plot the network
plot(plot_center_y ~ plot_center_x, coords, pch = 16, col = "blue")
for(i in 1:53){
  for(j in 1:53){
    if(i < j)
      lines(c(xy$x[i],xy$x[j]),c(xy$y[i],xy$y[j]),col = col.alpha("black",R_med[i, j]))
  }
}

m_2 <- map2stan(
  alist(
    y ~ dnorm(mu, sigma_res),
    mu = X * beta + gamma
  )
)



# dev incidence functions Hanski
dat <- data.frame(x = runif(3), y = runif(3), area = runif(3,1,10))
distance <- as.matrix(dist(dat[,1:2]))
area <- dat[,3]

isolation <- function(distance, area, alpha){
  n <- nrow(distance)
  isol <- rep(0,n)
  for(i in 1:n){
    isol[i] <- sum(area[-i] * exp(-alpha * distance[i,-i]))
  }
  return(isol)
}

saveRDS(isolation, "~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/scripts/isolation.r")

# get the needed data (distance between all forest patches and area of all forest patches)
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)

#set working directory
setwd("/home/lionel/Documents/PostDoc_Ghent/Data/GISdata/")

#now the raster as pixel sizes of 10 x 10m
fld <- raster(nrows=1200,ncols=2700,xmn=94700,xmx=121700,ymn=176000,ymx=188000,
              crs=CRS("+init=epsg:31370"))

#load land cover data
#lcc <- readOGR(dsn="LandCover/Shapefile/",layer="BwkHab")
#crop to the extent of treeweb
#lcc_treeweb <- crop(lcc,extent(94000,121584,175914,188764))
lcc_treeweb <- readOGR(dsn="LandCover/Shapefile/",layer="LandCover_TREEWEB")
#rename the land cover categories
forest <- c("^f","^q","^v")
#forest <- c("^sp","^sk","^sg","^sz","^sf","^so","^sm","^se","^ru","^vc","^va","^vf","^vn","^vm","^vo","^vt","^q","^f","^l","^p","^n","^kp","^kb","^kh")
toClass <- lcc_treeweb@data$EENH1

isForest <- unlist(lapply(forest,function(x) grep(x,toClass)))
lcc_treeweb@data$isForest <- 0
lcc_treeweb@data$isForest[isForest] <- 1
#only keep forest polygons
lcc_f <- subset(lcc_treeweb,isForest==1)
lcc_f2 <- spTransform(lcc_f,CRS("+init=epsg:31370"))
# save that
writeOGR(lcc_f,dsn = ".",layer = "forest_shp_daan2",driver = "ESRI Shapefile")

# re-read after Dgis manip
lcc_f3 <- readOGR("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/","treeweb_with_id_3")
lcc_f4 <- st_as_sf(lcc_f3)
lcc_f4$area <- st_area(lcc_f4)
lcc_f4[which(is.na(lcc_f4$fragment_i)),"fragment_i"] <- 93:(93+328)

lcc_f4 %>%
  group_by(fragment_i) %>%
  summarise(area = sum(area)) %>%
  st_cast() -> lcc_f5
# save this
write_sf(lcc_f5, dsn = "~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data",
         "treeweb_forests",driver="ESRI Shapefile",update=TRUE)

# save vector of areas
areas <- data.frame(fragment_id = lcc_f5$fragment_i,area=as.numeric(lcc_f5$area / 1e4))
write.table(areas,"~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/fragment_area.csv",row.names=FALSE)

# get distances
distt <- as.matrix(st_distance(lcc_f5))
distt <- distt / 1000
write.table(distt,"~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/fragment_dist.csv",row.names=FALSE)

nbb <- poly2nb(lcc_f)
summary(nbb)

lcc_f2 <- st_as_sf(lcc_f)
inters <- st_is_within_distance(lcc_f2,lcc_f2,dist=10)
lcc_u <- list()
mm <- as.matrix(inters)
# mm[lower.tri(mm)] <- FALSE

for(i in 1:nrow(mm)){
  if(sum(mm[i,]) == 0){
    next
  }
  else{
    id <- which(mm[i,])
    lcc_u[[length(lcc_u) + 1]] <- st_union(lcc_f2[id,])
    mm[,id] <- FALSE
  }
  print(paste("row:",i,"done"))
}

lcc_u2 <- st_multipolygon(lcc_u)

spP <- lapply(lcc_u, function(x) as(x,'Spatial'))

IDs <- sapply(spP, function(x)
  slot(slot(x, "polygons")[[1]], "ID"))
length(unique(IDs)) == length(spP)

Spol1 <- SpatialPolygons(lapply(1:length(spP), function(i) {
  Pol <- slot(spP[[i]], "polygons")[[1]]
  slot(Pol, "ID") <- as.character(i)
  Pol
}))

tmp <- st_as_sf(Spol1)
tmp$area <- st_area(tmp)

inters_red <- inters
feat_n <- 1:nrow(lcc_f2)



while(length(inters_red) != 0){
  to_union <- inters_red[[1]]
  feat_n <- feat_n[-which(feat_n %in % to_union)]
  inters_red <- inters_red[-to_union]
  lcc_u <- c(lcc_u,st_union(lcc_f2[to_union,]))
  print(paste("Poly:",to_union,"done"))
}

# import the plot coords
plots <- read.csv("~/Documents/PostDoc_Ghent/Data/sharepoint.ugent.be/general datafiles/Basisdataset_plots.csv")
pts <- SpatialPointsDataFrame(plots[,c("plot_center_x","plot_center_y")],data = data.frame(plot.id = plots$id_plot, fragment.id = plots$fragment_id),proj4string = CRS("+init=epsg:31370"))
writeOGR(pts,dsn = "~/Documents/PostDoc_Ghent/Data/GISdata/",layer = "plots_treeweb",driver = "ESRI Shapefile")

# isolation works now simulate some data for model testing
distt <- as.matrix(read.table("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/fragment_dist.csv",sep=" ",head=TRUE))
dimnames(distt) <- NULL
areas <- read.csv("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/fragment_area.csv",sep=" ")
isolation <- readRDS("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/scripts/isolation.r")


ii <- isolation(distt, areas$area, 0.05)
X <- runif(53,-2,2)
isol <- ii[plots$fragment_id]
isol <- scale(isol)
mu <- 1 + 2 * X + isol
y <- rnorm(53,mu,0.1)
X <- matrix(c(rep(1,53),X),ncol=2,byrow=FALSE)
frag_id <- plots$fragment_id

# create index data
id_m <- matrix(0,nrow=210,ncol=211,byrow=FALSE)
for(i in 1:211){
  id_m[,i] <- (1:211)[-i]
}

m <- stan(file = "~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/model/normal_isolation.stan",
          data = list(N=53,K=2,X=X,id_need=19,y=y,frag_id=frag_id,N_frag=420,
                      distt=distt,area=areas$area))

# use instead a proximity index
diag(distt) <- NA
distt[distt==0] <- NA
prox <- apply(distt, 1, function(x) sum(areas$area / (x ** 2),na.rm=TRUE))

X <- factor(sample(1:4,53,replace=TRUE))
proxi <- prox[plots$fragment_id]
proxi <- scale(proxi)
dat <- data.frame(X=X,proxi=proxi)
modmat <- model.matrix(~X+proxi,dat, contrasts = list(X = "contr.sum"))
X_n <- as.numeric(levels(X))[X]
mu <- modmat %*% c(1, 0.5, 1.3, -2, -1)
dat$y <- rnorm(53,mu,0.1)
#X <- matrix(c(rep(1,53),X,proxi[,1]),ncol=3,byrow=FALSE)
# dat <- data.frame(X = X, proxi = proxi, y = y)
# frag_id <- plots$fragment_id

# via brms
m <- brm(y ~ proxi + X, dat, family = gaussian())

# via stan
m <- stan("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/model/normal_proxi_brms.stan",
          data = list(n = 53, k = 5, n_mix = 4, y = dat$y, mix = X_n, proxi = dat$proxi, X = modmat))

# compute explained variation per elements
ss <- extract(m)
ss_mix <- ss$s_mix / (ss$s_mix + ss$s_proxi + ss$s_y)
ss_proxi <- ss$s_proxi / (ss$s_mix + ss$s_proxi + ss$s_y)
ss_y <- ss$s_y / (ss$s_mix + ss$s_proxi + ss$s_y)


# works!

# now turn to some real data with edge and composition effects
fragm <- read.csv("~/Documents/PostDoc_Ghent/Data/sharepoint.ugent.be/general datafiles/Fragmentation_plot.csv")
func <- read.csv("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/spatsynth_fun.csv")
## create the data
plots$rich <- ifelse(plots$spec_comp=="all",3,ifelse(nchar(as.character(plots$spec_comp))==4,1,2))
mix <- plots$spec_comp
mix_n <- as.numeric(plots$spec_comp)
plots$edge_eff <- scale(fragm$plot_edge100m)
plots$proxi <- scale(prox[plots$fragment_id])
# rich <- scale(plots$rich)
# frag_area <- scale(fragm$frag_area)

# some plotting
func %>%
  gather("fun","value",-id_plot) %>%
  left_join(plots[,c("id_plot","proxi","edge_eff","spec_comp")]) -> dd

ggplot(dd,aes(x=proxi,y=value))+
  geom_point(aes(color=spec_comp)) +
  stat_smooth(method="lm") +
  facet_wrap(~fun,scales="free_y")

dat <- data.frame(comp = mix, edge = edge_eff, proxi = proxi)
modmat <- model.matrix(~ comp + edge + proxi, dat, contrasts.arg = list(comp="contr.sum"))

# dat <- data.frame(comp = mix, proxi = proxi, edge = edge_eff, area = frag_area)
# modmat <- model.matrix(~ comp + proxi + edge + area, dat, contrasts.arg = list(comp="contr.sum"))


y <- func$Arth_div

m <- stan("normal_proxi_brms.stan",
          data = list(n = length(y), k = ncol(modmat), n_mix = length(unique(mix_n)), y = y, mix = mix_n, proxi = dat$proxi, edge = dat$edge, X = modmat))

m <- stan("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/model/normal_isolation.stan",
          data = list(N = length(y), K = ncol(modmat), id_need = 19, N_frag = 420,
                      X = modmat, y = y, frag_id = fragm$id_frag, distt = distt,
                      area = areas$area))
# compute explained variation per elements
ss <- extract(m)
ss_mix <- ss$s_mix / (ss$s_mix + ss$s_proxi + ss$s_edge + ss$s_y)
ss_proxi <- ss$s_proxi / (ss$s_mix + ss$s_proxi + ss$s_edge + ss$s_y)
ss_edge <- ss$s_edge / (ss$s_mix + ss$s_proxi + ss$s_edge + ss$s_y)
ss_y <- ss$s_y / (ss$s_mix + ss$s_proxi + ss$s_edge + ss$s_y)

# works ! 

# now go through all the functions and fit this model, removing NA rows
mm <- list()
setwd("Documents/PostDoc_Ghent/Spatialsynthesis_stuff/")
for(i in 2:22){
  y <- func[,i]
  y_nona <- y[!is.na(y)]
  if(length(y) != length(y_nona)){
    na <- which(!is.na(y))
    mix_nona <- mix_n[na]
    proxi_nona <- dat$proxi[na]
    edge_nona <- dat$edge[na]
    X_nona <- modmat[na,]
    #area_nona <- dat$rich[na]
    
    m <- stan("model/normal_proxi_brms.stan",
              data = list(n = length(y_nona), k = ncol(modmat), n_mix = length(unique(mix_n)), y = y_nona, mix = mix_nona, proxi = proxi_nona,
                          edge = edge_nona, X = X_nona))
    
  }
  else{
    m <- stan("model/normal_proxi_brms.stan",
              data = list(n = length(y), k = ncol(modmat), n_mix = length(unique(mix_n)), y = y, mix = mix_n, proxi = dat$proxi,
                          edge = dat$edge, X = modmat))
    
  }
  mm[names(func)[i]] <- m
  print(paste("function:",names(func)[i],"done"))
}

# grab and compute variance component
grab_ss <- function(model){
  ss <- rstan::extract(model)
  ss_mix <- ss$s_mix / (ss$s_mix + ss$s_proxi + ss$s_edge + ss$s_y)
  ss_proxi <- ss$s_proxi / (ss$s_mix + ss$s_proxi + ss$s_edge + ss$s_y)
  ss_edge <- ss$s_edge / (ss$s_mix + ss$s_proxi + ss$s_edge + ss$s_y)
  ss_y <- ss$s_y / (ss$s_mix + ss$s_proxi + ss$s_edge + ss$s_y)
  out <- data.frame(Variable = c("composition","proximity","edge","residuals"), SS = c(median(ss_mix),median(ss_proxi),median(ss_edge),median(ss_y)),
                    LCI = c(quantile(ss_mix,probs = 0.1),quantile(ss_proxi,probs = 0.1),quantile(ss_edge,probs=0.1), quantile(ss_y,probs=0.1)),
                    UCI = c(quantile(ss_mix,probs = 0.9),quantile(ss_proxi,probs = 0.9),quantile(ss_edge,probs=0.9), quantile(ss_y,probs=0.9)))
  return(out)
}

ss_dat <- ldply(mm, function(model) grab_ss(model))
ss_dat$.id <- factor(ss_dat$.id,
                     levels = c("C_stock","BS","pH","CN","P","Decomp",
                                "Biomass","Cover","Veg_div","GLI","LAI","P_germ","Seed_biom",
                                "Arth_div","Herbivory","Fit_spider","Spider_diet",
                                "Bird_smi","Breed_succ","Bird_div","Predation"))

ss_dat$X <- rep(1:21,each=4) + rep(c(-0.15,-0.07,0.07,0.15),times = 21)

ggplot(ss_dat,aes(x=.id,y=SS,fill=Variable)) +
  geom_bar(stat="identity") +
  coord_flip()


# another version
## get average per Variable
ss_dat %>%
  dplyr::group_by(Variable) %>%
  dplyr::summarise(SS = mean(SS)) -> ss_dd

gg_var <- ggplot(ss_dat,aes(x=X,y=SS,ymin=LCI,ymax=UCI,color=Variable)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  geom_hline(data=ss_dd,aes(yintercept=SS,color=Variable),linetype="dashed") +
  scale_x_continuous(breaks=1:21,labels = levels(ss_dat$.id)) +
  labs(x="Function",y="Proportion of variance explained")

ggsave("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/figures/variance_explained.png",gg_var)

# now look at edge effects
sds <- apply(func[,2:22],2,sd,na.rm=TRUE)

grab_edge <- function(index,model){
  post <- extract(mm[[index]])
  qq <- quantile(post$beta[,9] / sds[index],probs = c(0.1,0.5,0.9))
  out <- data.frame(Med = qq[2],LCI=qq[1],UCI=qq[3])
  return(out)
}

n <- names(mm)
names(n) <- names(mm)
edge_dat <- ldply(n,function(index,mm) grab_edge(index,mm))
edge_dat$.id <- factor(edge_dat$.id,
                     levels = c("C_stock","BS","pH","CN","P","Decomp",
                                "Biomass","Cover","Veg_div","GLI","LAI","P_germ","Seed_biom",
                                "Arth_div","Herbivory","Fit_spider","Spider_diet",
                                "Bird_smi","Breed_succ","Bird_div","Predation"))
gg_edge <- ggplot(edge_dat,aes(x=.id,y=Med,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  labs(x="Function",y="Edge effect (length of edge in a 100m radius)")

ggsave("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/figures/edge_slope.png",gg_edge)

# now look at proximity effects
grab_prox <- function(index,model){
  post <- extract(mm[[index]])
  qq <- quantile(post$beta[,8] / sds[index],probs = c(0.1,0.5,0.9))
  out <- data.frame(Med = qq[2],LCI=qq[1],UCI=qq[3])
  return(out)
}

n <- names(mm)
names(n) <- names(mm)
prox_dat <- ldply(n,function(index,mm) grab_prox(index,mm))
prox_dat$.id <- factor(prox_dat$.id,
                       levels = c("C_stock","BS","pH","CN","P","Decomp",
                                  "Biomass","Cover","Veg_div","GLI","LAI","P_germ","Seed_biom",
                                  "Arth_div","Herbivory","Fit_spider","Spider_diet",
                                  "Bird_smi","Breed_succ","Bird_div","Predation"))
gg_prox <- ggplot(prox_dat,aes(x=.id,y=Med,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  labs(x="Function",y="Proximity effect (positive values indicate function increase when closer to large forest fragments)")

ggsave("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/figures/prox_slope.png",gg_prox)

# now look at composition effect
grab_comp <- function(index){
  post <- extract(mm[[index]])
  qq <- apply(post$Beta_mix,2,function(x) quantile(x, probs = c(0.1,0.5,0.9)))
  out <- as.data.frame(t(qq))
  out$comp <- c("all","fsyl","qrob","qrob_fsyl","qrob_qrub","qrub","qrub_fsyl")
  out$X <- c(0.3,-0.3,-0.2,0,0.2,-0.1,0.1)
  names(out)[1:3] <- c("LCI","Med","UCI")
  return(out)
}

comp_dat <- ldply(n,function(index) grab_comp(index))

comp_dat$.id <- factor(comp_dat$.id,
                       levels = c("C_stock","BS","pH","CN","P","Decomp",
                                  "Biomass","Cover","Veg_div","GLI","LAI","P_germ","Seed_biom",
                                  "Arth_div","Herbivory","Fit_spider","Spider_diet",
                                  "Bird_smi","Breed_succ","Bird_div","Predation"))

comp_dat$X <- rep(1:21,each=7) + comp_dat$X
comp_dat$comp <- factor(comp_dat$comp,levels=c("fsyl","qrob","qrub","qrob_fsyl","qrub_fsyl","qrob_qrub","all"))  
comp_dat %>%
  group_by(.id) %>%
  summarise(Avg = mean(Med)) -> comp_dd
  
gg_comp <- ggplot(comp_dat,aes(x=comp,y=Med,ymin=LCI,ymax=UCI,color=comp)) +
  geom_linerange() +
  geom_point() +
  geom_hline(data=comp_dd,aes(yintercept=Avg),linetype="dashed") +
  facet_wrap(.~.id,scales = "free") +
  labs(x="Function",y="Composition effect") +
  theme(axis.text.x = element_text(angle=45,hjust=0.7))


# now try to add computation of explained variance
stan_c <- make_stancode(y ~ 0 + intercept + proxi + X, dat, family = gaussian(), prior = c(set_prior("normal(0,5)",class="b"),
                                                                            set_prior("cauchy(0,2)",class="sigma")),save_model = "normal_proxi_brms.stan")

lkj_corr <- function(l, eta){
  beta <- eta + (l-1)/2
  P <-  S <- matrix(0,nrow=l,ncol=l)
  diag(S) <- 1
  
  for(k in 1:(l-1)){
    beta <- beta - 0.5
    for(i in (k+1):l){
      P[k,i] <- rbeta(1, beta, beta)
      P[k,i] <- (P[k,i] - 0.5) * 2
      p <- P[k,i]
      if(k == 1){
        x <- 0
      }
      else{
        for(ll in seq((k-1),1,-1)){
          p = p * sqrt((1-P[ll,i] ** 2) * (1-P[ll,k]) ** 2) + P[ll,i] * P[ll,k]
        }
      }
      S[k,i] <- p
      S[i,k] <- p
    }
  }
  return(S)
}

# now go hierarchical by generating many variables
library(MASS)
N <- 100
K <- 3
J <- 10
L <- 2
jj <- sample(1:J, N, replace=TRUE)
X <- data.frame(int=1,x1=runif(N,-2,2),
                x2=runif(N,-2,2))
u <- data.frame(int=1,u1=runif(J,-2,2))

gamma <- matrix(c(1,0.5,-1,0.2,-0.3,0.5),ncol=K,byrow=TRUE)
u_gamma <- as.matrix(u) %*% gamma

omega <- lkj_corr(K, 2)

beta <- t(apply(u_gamma, 1, function(mu) mvrnorm(mu=mu,Sigma=omega)))

mu <- rep(NA,N)
for(n in 1:N){
  mu[n] <- as.matrix(X[n,]) %*% as.matrix(beta[jj[n],])
}

y <- rnorm(N, mu, 0.1)

m <- stan("model/normal_hier.stan",
          data = list(N=N,K=K,J=J,L=L,jj=jj,
                      x=X,u=u,y=y))


l <- 10 # number of groups
beta <- c(1, 2, -1) # the overall regression coeff
# generate correlation matrix between the varying effects
R <- lkj_corr(length(beta),1)
# the marginal standard deviations 
sd_beta <- c(0.1, 1, 2)
S <- sd_beta %*% t(sd_beta) * R
beta_gr <- mvrnorm(l, mu = beta, Sigma = S)

# the covariates
distt <- as.matrix(read.table("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/fragment_dist.csv",sep=";"))
dimnames(distt) <- NULL
areas <- read.csv("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/fragment_area.csv",sep=" ")
diag(distt) <- NA
distt[distt==0] <- NA
prox <- apply(distt, 1, function(x) sum(areas$area / x,na.rm=TRUE))

X <- runif(n,-2,2)
proxi <- prox[plots$fragment_id]
proxi <- scale(proxi)
dat <- data.frame(X = X, proxi = proxi)
modmat <- model.matrix(~ X + proxi, dat)

y_mat <- melt(apply(beta_gr, 1, function(B) modmat %*% B))
dat2 <- y_mat
dat2$X <- dat[y_mat$Var1,"X"]
dat2$proxi <- dat[y_mat$Var1,"proxi"]
names(dat2)[c(2,3)] <- c("Group","mu")
dat2$Group <- factor(paste0("G",dat2$Group))
dat2$y <- rnorm(530,mean = dat2$mu, sd = 0.1)

m <- brm(y ~ X + proxi + (X + proxi | Group), dat2, family = gaussian(), 
         prior = c(set_prior("normal(0,5)",class="b"),
                   set_prior("cauchy(0,2)",class="sigma")))


# try another approach with different spatial lags
## load back the vector data
vv <- read_sf("~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/data/treeweb_forests.shp")

#now the raster as pixel sizes of 10 x 10m
fld <- raster(nrows=1200,ncols=2700,xmn=94700,xmx=121700,ymn=176000,ymx=188000,
              crs=CRS("+init=epsg:31370"))

rr <- rasterize(vv,fld)
rr_v <- rr
values(rr_v)[!is.na(values(rr_v))] <- 1

# a function to get the surface area of forest within a certain buffer
get_forest_area <- function(raster, pts, buffer){
  ee <- extract(raster, pts, buffer = buffer)
  out <- ldply(ee, function(x) (sum(x, na.rm = TRUE) * prod(res(rr))) / 10000) # results in ha
  out$buffer <- buffer
  out$id <- 1:nrow(out)
  return(out)
}

dd <- ldply(seq(50,1000,20),function(x) get_forest_area(rr_v,pts,x))


spherical_decay <- function(x, phi){
  out <- ifelse(x > phi, 0,
                1-((3/2) * (x/phi)) + ((1/2) * (x/phi)**3))
  return(out)
}

extract(rr_v,plots[,c("plot_center_x","plot_center_y")],buffer=100)


