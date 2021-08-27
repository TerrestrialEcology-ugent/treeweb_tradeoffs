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
library(rstan)
rstan_options(auto_write = TRUE)
library(plyr)

N <- 100 # number of data points
K <- 3 # number of covariates
J <- 10 # number of groups
L <- 2 # number of group-level predictors
jj <- sample(1:J, N, replace=TRUE) # observation-level group index
X <- data.frame(int=1,x1=runif(N,-2,2),
                x2=runif(N,-2,2)) # obervation-level model matrix
u <- data.frame(int=1,u1=runif(J,-2,2)) # group-level model matrix

gamma <- matrix(c(1,0.5,-1,0.2,-0.3,0.5),ncol=K,byrow=TRUE) # group-level slopes
u_gamma <- as.matrix(u) %*% gamma # group-level betas

omega <- lkj_corr(K, 2) # correlation matrix between betas
diag(omega) <- c(1,2,3)

beta <- t(apply(u_gamma, 1, function(mu) mvrnorm(mu=mu,Sigma=omega))) # 

mu <- rep(NA,N)
for(n in 1:N){
  mu[n] <- as.matrix(X[n,]) %*% as.matrix(beta[jj[n],])
}

y <- rnorm(N, mu, 0.1)


m <- stan("model/normal_hier.stan",
          data = list(N=N,K=K,J=J,L=L,jj=jj,
                      x=X,u=u,y=y))

# explore model output
post <- extract(m)
## get group-level coefficient
gg <- adply(post$gamma, c(2,3), quantile, probs = c(0.1,0.5,0.9)) 
names(gg) <- c("GL_var", "Beta", "LCI", "Med", "UCI")
gg$GL_var <- rep(c("Avg","u1"),3)
gg$Beta <- rep(c("Int","x1","x2"),each=2)
gg$X <- rep(1:3,each=2) + rep(c(-0.1,0.1),3)

ggplot(gg,aes(x=X,y=Med,ymin=LCI,ymax=UCI,color=GL_var)) +
  geom_linerange() +
  geom_point() +
  scale_x_continuous(breaks = 1:3,labels = c("Int","x1","x2"))

## now get the betas
bb <- adply(post$beta, c(2,3), quantile, probs = c(0.1,0.5,0.9)) 
names(bb) <- c("groups","predictor","LCI","Med","UCI")
bb$groups <- paste0("Fun",bb$groups)
bb$predictor <- rep(c("Int","x1","x2"),each=10)

ggplot(bb,aes(x=groups,y=Med,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  facet_grid(~predictor) +
  coord_flip()

## now get the finite sd
ss <- adply(cbind(post$sd_y,post$sd_x),2,quantile,probs=c(0.1,0.5,0.9))

# fit this model to real data
## create the response vector
func %>%
  select(id_plot:Bird_div) %>%
  gather("Function","value",-id_plot) %>%
  group_by(Function) %>%
  mutate(value_std = scale(value))-> y_real

y_real[is.na(y_real$value),"value_std"] <- 0

u_real <- data.frame(Function = unique(y_real$Function),
                     Compartment = c(rep("Soil",6),rep("Veg",7),rep("Arth",2),"Bird","Arth","Arth",rep("Bird",3)))

uu_real <- model.matrix(~Compartment, u_real, contrasts.arg = list(Compartment = "contr.sum"))

x_real <- model.matrix(~spec_comp + edge_eff + proxi, plots, contrasts.arg = list(spec_comp = "contr.sum"))
xx_real <- x_real[rep(1:53,21),]

jj <- rep(1:21,each=53)
N <- length(jj)
K <- ncol(xx_real)
J <- nrow(uu_real)
L <- ncol(uu_real)

m <- stan("model/normal_hier.stan",
          data = list(N=N,K=K,J=J,L=L,jj=jj,x=xx_real,u=uu_real,y=y_real$value_std),
          control = list(adapt_delta = 0.99))

## get group-level coefficient
gg <- adply(post$gamma, c(2,3), quantile, probs = c(0.1,0.5,0.9)) 
names(gg) <- c("GL_var", "Beta", "LCI", "Med", "UCI")
gg$GL_var <- rep(c("Avg","u1"),3)
gg$Beta <- rep(c("Int","x1","x2"),each=2)
gg$X <- rep(1:3,each=2) + rep(c(-0.1,0.1),3)

ggplot(gg,aes(x=X,y=Med,ymin=LCI,ymax=UCI,color=GL_var)) +
  geom_linerange() +
  geom_point() +
  scale_x_continuous(breaks = 1:3,labels = c("Int","x1","x2"))

## now get the betas
bb <- adply(post$beta, c(2,3), quantile, probs = c(0.1,0.5,0.9)) 
names(bb) <- c("groups","predictor","LCI","Med","UCI")
bb$groups <- paste0("Fun",bb$groups)
bb$predictor <- rep(c("Int","x1","x2"),each=10)

ggplot(bb,aes(x=groups,y=Med,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  facet_grid(~predictor) +
  coord_flip()

## now get the finite sd
ss <- adply(cbind(post$sd_y,post$sd_x),2,quantile,probs=c(0.1,0.5,0.9))

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

############ start brms analysis #########
library(plyr)
library(tidyverse)
library(brms)
library(circlize)

# try brms hierarchical stuff
func <- read.csv("data/spatsynth_fun.csv")
div_all <- read.csv("data/community_data_formatted.csv")
fragm <- read.csv("data/Fragmentation_plot.csv")
plots <- read.csv("data/spatsynth_expl.csv")


pred_dat <- cbind(fragm[,c("id_plot","frag_prox","plot_edge100m")],plots[,c("speccomb", "total_all_ba")])
names(pred_dat)[4] <- "speccomb"

# create a dataset for the function data
func %>%
  select(id_plot,C_stock,Decomp,Biomass,P_germ,Herbivory,Predation,Bird_smi,frass) %>%
  gather("Function","value",-id_plot) %>%
  group_by(Function) %>%
  mutate(value_std = scale(value)) %>%
  mutate(value_std = ifelse(is.na(value_std), 0, value_std)) %>%
  left_join(pred_dat, by = "id_plot") %>%
  mutate(prox_std = scale(frag_prox), edge_std = scale(plot_edge100m)) -> brm_fun


m_func <- brm(value_std ~ 0 + intercept + speccomb + prox_std + edge_std + (0 + intercept + speccomb + prox_std + edge_std | Function), 
         data = brm_fun, family = gaussian,
         prior = set_prior("normal(0,2)",class="b"))

# some checks
pp_check(m_func) # good
hist(rhat(m_func)) # good
hist(neff_ratio(m_func)) # ok

# posterior draws
pp <- posterior_samples(m_func)
ff <- as.data.frame(fixef(m_func)) # population-level effects
ff$Eff <- c("all","fsyl","fsyl_qrob","fsyl_qrub","qrob","qrob_qrub","qrub","proximity","edge")
rr <- ranef(m_func)
func_n <- unique(brm_dat_func$Function)[order(unique(brm_dat_func$Function))]

## look at population-level and group-level effects
# first fragment effect
ee <- pp$b_edge_std + pp[,138:145]
ee %>%
  gather("Fun","value") %>%
  group_by(Fun) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(Fun = func_n, Eff = "edge") -> ee_dd

prox <- pp$b_prox_std + pp[,130:137]
prox %>%
  gather("Fun","value") %>%
  group_by(Fun) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(Fun = func_n, Eff = "proximity") -> prox_dd

frag_eff <- rbind(ee_dd,prox_dd)

gg_frag_fun <- ggplot(frag_eff,aes(x=Fun,y=Estimate,ymin=Q2.5,ymax=Q97.5)) +
  geom_linerange() +
  geom_point() +
  geom_hline(data=ff[8:9,],aes(yintercept=Estimate)) +
  geom_hline(data=ff[8:9,],aes(yintercept=Q2.5),linetype="dashed") +
  geom_hline(data=ff[8:9,],aes(yintercept=Q97.5),linetype="dashed") +
  coord_flip() +
  facet_grid(~Eff)

ggsave("figures/frag_fun.png",gg_frag_fun)

# now composition effect
# poulation-level effect
ff_comp <- cbind(pp$b_intercept,pp$b_intercept + pp[,2:7])
ff_comp %>%
  gather("comp","value") %>%
  group_by(comp) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(speccomb = c("fsyl","fsyl_qrob","fsyl_qrub","qrob","qrob_qrub","qrub","all")) -> ff_comp_dd

# group-level effects
newdat <- expand.grid(speccomb = unique(brm_dat_func$speccomb), edge_std = 0, prox_std = 0, Function = func_n)

out_comp <- cbind(newdat,fitted(m_func, newdata = newdat))

out <- NULL
for(i in 1:8){
  int <- 66:73
  ids <- seq(int[i],129,8)
  tmp <- pp[,ids[1]] + pp[,ids[2:8]] + pp[,1] + pp[,2:8]
  tmp %>%
    gather("comp","value") %>%
    group_by(comp) %>%
    summarise(Estimate = median(value),
              Q2.5 = quantile(value, probs = 0.025),
              Q97.5 = quantile(value, probs = 0.975)) %>%
    mutate(comp = c("all","fsyl_qrob","fsyl_qrub","fsyl","qrob_qrub","qrob","qrub")) %>%
    mutate(func = func_n[i]) -> tmp_dd
  out <- rbind(out,tmp_dd)
}

comp_eff <- ggplot(out_comp,aes(x=Function,y=Estimate,ymin=Q2.5,ymax=Q97.5)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  facet_wrap(~speccomb) +
  geom_hline(data=ff_comp_dd, aes(yintercept = Estimate)) +
  geom_hline(data=ff_comp_dd, aes(yintercept = Q2.5),linetype="dashed") +
  geom_hline(data=ff_comp_dd, aes(yintercept = Q97.5),linetype="dashed") 

ggsave("figures/comp_fun.png",comp_eff)

# compute finite-population variance
newdat <- expand.grid(Function = func_n, speccomb = unique(brm_dat_func$speccomb),edge_std=0,prox_std=0)
pp_comp <- posterior_linpred(m_func,newdata=newdat)
sd_comp <- as.data.frame(t(apply(sapply(seq(1,105,15),function(x) apply(pp_comp[,x:(x+14)],1,sd)),2,quantile,probs=c(0.025,0.5,0.975)))) # palme d'or code le plus incomprehensible 2019
sd_comp$Var <-  unique(brm_dat_func$speccomb)
names(sd_comp)[1:3] <- c("Q2.5", "Estimate","Q97.5")

# out_sd <- data.frame(All = pp$b_intercept + pp[,56])
# 
# 
# 
# for(i in 1:14){
#   int <- 56:69
#   ids <- seq(int[i],153,14)
#   tmp <- cbind(pp[,ids[1]] + pp[,1],pp[,ids[1]] + pp[,ids[2:7]] + pp[,1] + pp[,2:7])
#   tmp %>%
#     gather("comp","value") %>%
#     group_by(comp) %>%
#     summarise(Estimate = median(value),
#               Q2.5 = quantile(value, probs = 0.025),
#               Q97.5 = quantile(value, probs = 0.975)) %>%
#     mutate(comp = c("all","fsyl_qrob","fsyl_qrub","fsyl","qrob_qrub","qrob","qrub")) %>%
#     mutate(func = func[i]) -> tmp_dd
#   out <- rbind(out,tmp_dd)
# }
# 
# comp <- list(all = post[,56:76])
# ff <- seq(77,202,21)
# comp_n <- c("fsyl","qrob","qrob_fsyl","qrob_qrub","qrub","qrub_fsyl")
# for(i in 1:6){
#     comp[[comp_n[i]]] <- comp$all + post[,ff[i]:(ff[i] + 20)]
# }
# 
# comp_dd <- ldply(comp)
# comp_dd$iter <- rep(1:4000,7)
# comp_dd %>%
#   gather("Function","value",-.id, -iter) -> comp_dd2
# 
# comp_dd2 %>%
#   group_by(.id, iter) %>%
#   summarise(SD = sd(value)) -> sd_comp
# 
# sd_comp %>%
#   group_by(.id) %>%
#   summarise(LCI=quantile(SD,0.1),Med=median(SD),UCI=quantile(SD,0.9)) %>%
#   rename(Var=.id) -> sd_all
# 
# ggplot(sd_comp,aes(x=SD,fill=.id)) +
#   geom_density(alpha=0.2)

# post <- posterior_samples(m_func)
ee <- pp[,201:215]
ee_sd <- apply(ee,1,sd)
qq_ee <- quantile(ee_sd,probs=c(0.025,0.5,0.975))
ee_dd <- data.frame(Q2.5=qq_ee[1],Estimate=qq_ee[2],Q97.5=qq_ee[3],Var="edge")

prox <- pp[,186:200]
pp_sd <- apply(prox,1,sd)
qq_pp <- quantile(pp_sd,probs=c(0.025,0.5,0.975))
pp_dd <- data.frame(Q2.5=qq_pp[1],Estimate=qq_pp[2],Q97.5=qq_pp[3],Var="prox")


# get residuals variation
fitt <- fitted(m_func,summary = FALSE)
resid <- apply(fitt,1, function(y_hat) brm_dat_func$value_std - y_hat)
qq_resid <- quantile(apply(resid,1,sd),probs=c(0.025,0.5,0.975))
sd_resid <- data.frame(Q2.5=qq_pp[1],Estimate=qq_pp[2],Q97.5=qq_pp[3],Var="resid")

# figure
sd_all <- rbind(sd_comp,ee_dd,pp_dd,sd_resid)

gg_sd <- ggplot(sd_all,aes(x=Var,y=Estimate,ymin=Q2.5,ymax=Q97.5)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  labs(y = "Finite population deviation (with 95% CrI)",x="Predictors")


ggsave("figures/finite_sd.png",gg_sd)

## now for diversity
div_all %>%
  select(id_plot,Sha_bat,Sha_bird,Sha_spider,Sha_herb,Sha_veg,Sha_cara,Sha_iso,Sha_dip) %>%
  mutate(id_plot = as.numeric(id_plot)) %>%
  gather("group","diversity",-id_plot) %>%
  mutate(group = gsub("Sha_","",group)) %>%
  left_join(pred_dat, by = "id_plot") %>%
  mutate(prox_std = scale(frag_prox), edge_std = scale(plot_edge100m)) -> brm_div

# a big data frame
div_all %>%
  mutate(id_plot = as.numeric(id_plot)) %>%
  gather("group","diversity",-id_plot) %>%
  separate(group,c("index","group")) %>% 
  left_join(pred_dat, by = "id_plot") %>%
  mutate(prox_std = scale(frag_prox), edge_std = scale(plot_edge100m)) -> div_dd

# some plotting
ggplot(div_dd,aes(x=prox_std,y=diversity)) +
  geom_point() +
  stat_smooth(method="lm") +
  facet_wrap(index~group,scales = "free_y")

ggplot(div_dd,aes(x=edge_std,y=diversity)) +
  geom_point() +
  stat_smooth(method="lm") +
  facet_wrap(index~group,scales = "free_y")

ggplot(div_dd,aes(x=speccomb,y=diversity)) +
  geom_boxplot() +
  facet_wrap(index~group,scales="free_y")

m_div <- brm(diversity ~ 0 + intercept + speccomb + prox_std + edge_std + (0 + intercept + speccomb + prox_std + edge_std | group), 
             brm_div, prior = set_prior("normal(0,2)",class = "b"), family = gaussian, control = list(adapt_delta = 0.99))


# some checks
pp_check(m_div) # not so optimal some positive bias
hist(rhat(m_div)) # good
hist(neff_ratio(m_div)) # some parameter have low n_eff ratio

# posterior draws
pp <- posterior_samples(m_div)
ff <- as.data.frame(fixef(m_div)) # population-level effects
ff$Eff <- c("all","fsyl","fsyl_qrob","fsyl_qrub","qrob","qrob_qrub","qrub","proximity","edge")
rr <- ranef(m_func)
grp_n <- unique(brm_div$group)[order(unique(brm_div$group))]

## look at population-level and group-level effects
# first fragment effect
ee <- pp$b_edge_std + pp[,138:145]
ee %>%
  gather("group","value") %>%
  group_by(group) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(group = grp_n, Eff = "edge") -> ee_dd

prox <- pp$b_prox_std + pp[,130:137]
prox %>%
  gather("group","value") %>%
  group_by(group) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(group = grp_n, Eff = "proximity") -> prox_dd

frag_eff <- rbind(ee_dd,prox_dd)

gg_frag_div <- ggplot(frag_eff,aes(x=group,y=Estimate,ymin=Q2.5,ymax=Q97.5)) +
  geom_linerange() +
  geom_point() +
  geom_hline(data=ff[8:9,],aes(yintercept=Estimate)) +
  geom_hline(data=ff[8:9,],aes(yintercept=Q2.5),linetype="dashed") +
  geom_hline(data=ff[8:9,],aes(yintercept=Q97.5),linetype="dashed") +
  coord_flip() +
  facet_grid(~Eff)

ggsave("figures/frag_div.png",gg_frag_div)

# now composition effect
# poulation-level effect
ff_comp <- cbind(pp$b_intercept,pp$b_intercept + pp[,2:7])
ff_comp %>%
  gather("comp","value") %>%
  group_by(comp) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(speccomb = c("fsyl","fsyl_qrob","fsyl_qrub","qrob","qrob_qrub","qrub","all")) -> ff_comp_dd

# group-level effects
newdat <- expand.grid(speccomb = unique(brm_div$speccomb), edge_std = 0, prox_std = 0, group = grp_n)

out_comp <- cbind(newdat,fitted(m_div, newdata = newdat))

out <- NULL
for(i in 1:8){
  int <- 66:73
  ids <- seq(int[i],129,8)
  tmp <- pp[,ids[1]] + pp[,ids[2:8]] + pp[,1] + pp[,2:8]
  tmp %>%
    gather("comp","value") %>%
    group_by(comp) %>%
    summarise(Estimate = median(value),
              Q2.5 = quantile(value, probs = 0.025),
              Q97.5 = quantile(value, probs = 0.975)) %>%
    mutate(comp = c("all","fsyl_qrob","fsyl_qrub","fsyl","qrob_qrub","qrob","qrub")) %>%
    mutate(group = grp_n[i]) -> tmp_dd
  out <- rbind(out,tmp_dd)
}

comp_eff <- ggplot(out_comp,aes(x=group,y=Estimate,ymin=Q2.5,ymax=Q97.5)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  facet_wrap(~speccomb) +
  geom_hline(data=ff_comp_dd, aes(yintercept = Estimate)) +
  geom_hline(data=ff_comp_dd, aes(yintercept = Q2.5),linetype="dashed") +
  geom_hline(data=ff_comp_dd, aes(yintercept = Q97.5),linetype="dashed") 

ggsave("figures/comp_div.png",comp_eff)



post_div <- posterior_samples(m_div)
# plot edge effect per group
edge_eff <- as.data.frame(t(apply(post_div$b_edge_std + post_div[,96:100],2,quantile,probs=c(0.1,0.5,0.9))))
edge_eff$group <- c("bat","bird","herb","spider","veg")
names(edge_eff)[1:3] <- c("LCI","Med","UCI")

ggplot(edge_eff,aes(x=group,y=Med,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  geom_hline(yintercept =  0, linetype = "dashed")

# plot prox effect per group
prox_eff <- as.data.frame(t(apply(post_div$b_prox_std + post_div[,91:95],2,quantile,probs=c(0.1,0.5,0.9))))
prox_eff$group <- c("bat","bird","herb","spider","veg")
names(prox_eff)[1:3] <- c("LCI","Med","UCI")

ggplot(prox_eff,aes(x=group,y=Med,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  geom_hline(yintercept =  0, linetype = "dashed")

newdat <- expand.grid(prox_std = seq(-1.4,3.5,length=10),edge_std=0,speccomb = "all", group = c("bird","bat","spider","herb","veg"))
newdat$fitted <- fitted(m_div,newdata=newdat)[,1]

newdat %>%
  group_by(prox_std) %>%
  dplyr::summarise(fit = sum(fitted * c(10, 8, 5, 1, 8)) / sum(c(10,8,5,1,8))) -> dd



ggplot(dd,aes(x=prox_std,y=fit)) +
  geom_line()


newdat <- expand.grid(prox_std = 0,edge_std=seq(-1.15,2,length=10),speccomb = "all", group = c("bird","bat","spider","herb","veg"))
newdat$fitted <- fitted(m_div,newdata=newdat)[,1]

newdat %>%
  group_by(edge_std) %>%
  dplyr::summarise(fit = sum(fitted * c(10, 8, 5, 1, 8)) / sum(c(10,8,5,1,8))) -> dd

ggplot(dd,aes(x=prox_std,y=fit)) +
  geom_line()

newdat <- expand.grid(prox_std = 0, edge_std = 0, speccomb = c("all","fsyl","qrob","qrub","fsyl_qrob","fsyl_qrub","qrob_qrub"),group = c("bird","bat","spider","herb","veg","iso","dip","cara"))
newdat <- cbind(newdat,fitted(m_div,newdata=newdat))

newdat %>%
  group_by(speccomb) %>%
  dplyr::summarise(fit = sum(Estimate * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10,8,5,1,8,5,5,5)),
                   LCI = sum(Q2.5 * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10, 8, 5, 1, 8, 5, 5, 5)),
                   UCI = sum(Q97.5 * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10, 8, 5, 1, 8, 5, 5, 5))) -> comp_div

# make an ordered factor
comp_div$speccomb <- factor(comp_div$speccomb,levels = arrange(comp_div,fit)$speccomb)
# plot
p_comp_div <- ggplot(comp_div,aes(x=speccomb,y=fit,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  labs(x="Tree species composition (sorted)",y="Multidiversity index")

# compare with functions
newdat <- expand.grid(prox_std = 0, edge_std = 0, speccomb = c("all","fsyl","qrob","qrub","fsyl_qrob","fsyl_qrub","qrob_qrub"),
                      Function = unique(brm_dat_func$Function))
newdat <- cbind(newdat,fitted(m_func,newdata=newdat))
impA <- read.csv("../Synthesis_stuff/treeweb_synthesis/data/synthesis_importance_scores.csv", sep =" ")$Importance_weight
impA <- c(impA[c(1:7,10:13,15,17,20)],5)
newdat %>%
  group_by(speccomb) %>%
  dplyr::summarise(fit = sum(Estimate * impA) / sum(impA),
                   LCI = sum(Q2.5 * impA) / sum(impA),
                   UCI = sum(Q97.5 * impA) / sum(impA)) -> comp_func

comp_func$speccomb <- factor(comp_func$speccomb,levels = arrange(comp_func,fit)$speccomb)
# plot
p_comp_fun <- ggplot(comp_func,aes(x=speccomb,y=fit,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  labs(x="Tree species composition (sorted)",y="Multifunctionality index")

# now do the same for proximity
newdat <- expand.grid(prox_std = seq(-1.4,3.5, length=10),edge_std = 0, speccomb = "all",group = unique(brm_div$group))
newdat <- cbind(newdat,fitted(m_div,newdata = newdat))


newdat %>%
  group_by(prox_std) %>%
  dplyr::summarise(fit = sum(Estimate * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10,8,5,1,8,5,5,5)),
                   LCI = sum(Q2.5 * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10, 8, 5, 1, 8, 5, 5, 5)),
                   UCI = sum(Q97.5 * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10, 8, 5, 1, 8, 5, 5, 5))) -> prox_div

# plot
p_prox_div <- ggplot(prox_div,aes(x=prox_std,y=fit,ymin=LCI,ymax=UCI)) +
  geom_ribbon(alpha=0.1) +
  geom_line()  +
  labs(x="Proximity index (from isolated to well-connected)",y="Multidiversity index")

# for the functions
newdat <- expand.grid(prox_std = seq(-1.4,3.5, length=10),edge_std = 0, speccomb = "all",Function = unique(brm_dat_func$Function))
newdat <- cbind(newdat,fitted(m_func,newdata = newdat))


newdat %>%
  group_by(prox_std) %>%
  dplyr::summarise(fit = sum(Estimate * impA) / sum(impA),
                   LCI = sum(Q2.5 * impA) / sum(impA),
                   UCI = sum(Q97.5 * impA) / sum(impA)) -> prox_fun

# plot
p_prox_fun <- ggplot(prox_fun,aes(x=prox_std,y=fit,ymin=LCI,ymax=UCI)) +
  geom_ribbon(alpha=0.1) +
  geom_line()  +
  labs(x="Proximity index (from isolated to well-connected)",y="Multifunctionality index")

# now for edge
newdat <- expand.grid(prox_std = 0,edge_std = seq(-1.13,2,length=10), speccomb = "all",group = unique(brm_div$group))
newdat <- cbind(newdat,fitted(m_div,newdata = newdat))


newdat %>%
  group_by(edge_std) %>%
  dplyr::summarise(fit = sum(Estimate * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10,8,5,1,8,5,5,5)),
                   LCI = sum(Q2.5 * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10, 8, 5, 1, 8, 5, 5, 5)),
                   UCI = sum(Q97.5 * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10, 8, 5, 1, 8, 5, 5, 5))) -> edge_div

# plot
p_edge_div <- ggplot(edge_div,aes(x=edge_std,y=fit,ymin=LCI,ymax=UCI)) +
  geom_ribbon(alpha=0.1) +
  geom_line()  +
  labs(x="Edge amount in a buffer in 100m",y="Multidiversity index")

# for the functions
newdat <- expand.grid(prox_std = 0,edge_std = seq(-1.13,2,length=10), speccomb = "all",Function = unique(brm_dat_func$Function))
newdat <- cbind(newdat,fitted(m_func,newdata = newdat))


newdat %>%
  group_by(edge_std) %>%
  dplyr::summarise(fit = sum(Estimate * impA) / sum(impA),
                   LCI = sum(Q2.5 * impA) / sum(impA),
                   UCI = sum(Q97.5 * impA) / sum(impA)) -> edge_fun

# plot
p_edge_fun <- ggplot(edge_fun,aes(x=edge_std,y=fit,ymin=LCI,ymax=UCI)) +
  geom_ribbon(alpha=0.1) +
  geom_line()  +
  labs(x="Edge amount in a buffer of 100m",y="Multifunctionality index")

# all together
grid.arrange(p_comp_fun,p_comp_div,p_prox_fun,p_prox_div,p_edge_fun,p_edge_div,ncol=2)

# try to find optimum conditions
speccomb <- c("all","fsyl","qrob","qrub","fsyl_qrob","fsyl_qrub","qrob_qrub")
likfun_div <- function(par,comp="all"){
  edge <- as.numeric(par[1])
  prox <- as.numeric(par[2])
  
  newdat_div <- expand.grid(prox_std = prox, edge_std = edge, speccomb = comp, group = unique(brm_div$group))
  newdat_fun <- expand.grid(prox_std = prox, edge_std = edge, speccomb = comp, Function = unique(brm_dat_func$Function))
  
  p_div <- fitted(m_div,newdata=newdat_div)
  p_fun <- fitted(m_func, newdata = newdat_fun)
  
  mdiv <- sum(p_div[,1] * c(10, 8, 5, 1, 8, 5, 5, 5)) / sum(c(10, 8, 5, 1, 8, 5, 5, 5)) #multidiv
  mfun <- sum(p_fun[,1] * impA) / sum(impA)
  
  return(c(mdiv,mfun,mdiv + mfun))
}

trying <- expand.grid(edge = seq(-1.2,2,length=20),prox = seq(-1.4,3.5,length=20),comp = c("fsyl","qrob","qrub","fsyl_qrob","fsyl_qrub","qrob_qrub","all"))
trying$mdiv <- 0
trying$mfun <- 0
trying$mall <- 0

for(i in 1:nrow(trying)){
  tmp <- likfun_div(trying[i,1:2],trying[i,3])
  trying$mdiv[i] <- tmp[1]
  trying$mfun[i] <- tmp[2]
  trying$mall[i] <- tmp[3]
  print(i)
}

ggplot(trying,aes(x=edge,y=prox,fill=mdiv)) +
  geom_tile() +
  facet_wrap(~comp) +
  scale_fill_viridis()


ggplot(trying,aes(x=edge,y=prox,fill=mfun)) +
  geom_tile() +
  facet_wrap(~comp) +
  scale_fill_viridis()

ggplot(trying,aes(x=edge,y=prox,fill=mall)) +
  geom_tile() +
  facet_wrap(~comp) +
  scale_fill_viridis()


## now for abundance
div_all %>%
  select(id_plot,Abun_bat,Abun_bird,Abun_spider,Abun_herb,Abun_veg,Abun_cara,Abun_iso,Abun_dip) %>%
  mutate(Abun_veg = round(Abun_veg,0)) %>%
  #mutate_at(2:9,scale) %>%
  mutate(id_plot = as.numeric(id_plot)) %>%
  gather("group","abundance",-id_plot) %>%
  mutate(group = gsub("Abun_","",group)) %>%
  left_join(pred_dat, by = "id_plot") %>%
  mutate(prox_std = scale(frag_prox), edge_std = scale(plot_edge100m)) -> brm_abu

# the model
m_abu <- brm(abundance ~ speccomb + prox_std + edge_std + (speccomb + prox_std + edge_std | group),
             brm_abu, family = negbinomial, prior = set_prior("normal(0,2)",class="b"),control = list(adapt_delta=0.99))

# some checks
pp_check(m_abu) # looks ok
hist(rhat(m_div)) # good
hist(neff_ratio(m_div)) # some parameter have low n_eff ratio

# posterior draws
pp <- posterior_samples(m_abu)
ff <- as.data.frame(fixef(m_abu)) # population-level effects
ff$Eff <- c("all","fsyl","fsyl_qrob","fsyl_qrub","qrob","qrob_qrub","qrub","proximity","edge")

grp_n <- unique(brm_div$group)[order(unique(brm_div$group))]

## look at population-level and group-level effects
# first fragment effect
ee <- pp$b_edge_std + pp[,120:127]
ee %>%
  gather("group","value") %>%
  group_by(group) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(group = grp_n, Eff = "edge") -> ee_dd

prox <- pp$b_prox_std + pp[,112:119]
prox %>%
  gather("group","value") %>%
  group_by(group) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(group = grp_n, Eff = "proximity") -> prox_dd

frag_eff <- rbind(ee_dd,prox_dd)

gg_frag_abu <- ggplot(frag_eff,aes(x=group,y=Estimate,ymin=Q2.5,ymax=Q97.5)) +
  geom_linerange() +
  geom_point() +
  geom_hline(data=ff[8:9,],aes(yintercept=Estimate)) +
  geom_hline(data=ff[8:9,],aes(yintercept=Q2.5),linetype="dashed") +
  geom_hline(data=ff[8:9,],aes(yintercept=Q97.5),linetype="dashed") +
  coord_flip() +
  facet_grid(~Eff)

ggsave("figures/frag_abu.png",gg_frag_abu)

# now composition effect
# poulation-level effect
ff_comp <- cbind(pp$b_Intercept,pp$b_Intercept + pp[,2:7])
ff_comp %>%
  gather("comp","value") %>%
  group_by(comp) %>%
  summarise(Estimate = median(value),
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975)) %>%
  mutate(speccomb = c("fsyl","fsyl_qrob","fsyl_qrub","qrob","qrob_qrub","qrub","all")) -> ff_comp_dd

# group-level effects
newdat <- expand.grid(speccomb = unique(brm_div$speccomb), edge_std = 0, prox_std = 0, group = grp_n)

out_comp <- cbind(newdat,fitted(m_div, newdata = newdat))

out <- NULL
for(i in 1:8){
  int <- 56:63
  ids <- seq(int[i],111,8)
  tmp <- cbind(pp[,1] + pp[,ids[1]],pp[,ids[1]] + pp[,ids[2:7]] + pp[,1] + pp[,2:7])
  tmp %>%
    gather("comp","value") %>%
    group_by(comp) %>%
    summarise(Estimate = median(value),
              Q2.5 = quantile(value, probs = 0.025),
              Q97.5 = quantile(value, probs = 0.975)) %>%
    mutate(comp = c("all","fsyl_qrob","fsyl_qrub","fsyl","qrob_qrub","qrob","qrub")) %>%
    mutate(group = grp_n[i]) -> tmp_dd
  out <- rbind(out,tmp_dd)
}

comp_eff <- ggplot(out_comp,aes(x=group,y=Estimate,ymin=Q2.5,ymax=Q97.5)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  facet_wrap(~speccomb) +
  geom_hline(data=ff_comp_dd, aes(yintercept = Estimate)) +
  geom_hline(data=ff_comp_dd, aes(yintercept = Q2.5),linetype="dashed") +
  geom_hline(data=ff_comp_dd, aes(yintercept = Q97.5),linetype="dashed") 

ggsave("figures/comp_abu.png",comp_eff)


# some mediation, pcsem-type analysis
## put together a data frame
func %>%
  dplyr::select(id_plot,C_stock,Biomass,frass,Bird_smi,Decomp,P_germ,Herbivory,Predation) %>%
  gather(fun, value,-id_plot) %>%
  group_by(fun) %>%
  mutate(value_std = scale(value)) %>%
  ungroup() %>%
  left_join(pred_dat,by="id_plot") %>%
  mutate(edge = scale(plot_edge100m),prox=scale(frag_prox),rich=scale(rich))-> fun_dd

div_all %>%
  dplyr::select(id_plot,Sha_veg,Sha_herb,Sha_iso,Sha_dip,Sha_cara,Sha_spider,Sha_bird,Sha_bat) %>%
  gather(grp, diversity,-id_plot) %>%
  group_by(grp) %>%
  mutate(div_std = scale(diversity)) %>%
  separate(grp,c("drop","grp")) %>%
  ungroup() -> div_dd

func %>%
  dplyr::select(id_plot,C_stock,Biomass,frass,Bird_smi,Decomp,P_germ,Herbivory,Predation) %>%
  gather(fun, value,-id_plot) %>%
  group_by(fun) %>%
  mutate(value_std = scale(value)) %>%
  ungroup() %>%
  spread(fun,value)

fun_v$Sha_veg <- div_all$Sha_veg

bf_veg <- bf(Sha_veg ~ rich + edge + prox)
bf_treebiom <- bf(Biomass ~ rich + edge + prox + Sha_veg)

bform <- bf_veg + bf_treebiom + set_rescor(FALSE)
fit <- brm(bform,fun_v)

mm <- brm(mvbind(Biomass,Bird_smi,C_stock,Decomp,frass,Herbivory,P_germ,Predation) ~ rich + edge + prox, data = fun_v)
  

# build some models separate for each div components
expl <- data.frame(id_plot = fragm$id_plot, edge = fragm$plot_edge100m, prox = fragm$frag_prox, speccomp = plots$spec_comp)

div_all %>%
  gather("variable","value", -id_plot) %>%
  separate(variable,c("metric","taxa")) %>%
  spread(metric,value) %>%
  left_join(expl, by = "id_plot") -> div_expl

# effect on abun
ggplot(div_expl, aes(x = edge, y = Abun)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~taxa, scales = "free") +
  scale_y_log10()

ggplot(div_expl, aes(x = prox, y = Abun)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~taxa, scales = "free") +
  scale_y_log10()

ggplot(div_expl, aes(x = speccomp, y = Abun)) +
  geom_boxplot() +
  facet_wrap(~taxa, scales = "free") +
  scale_y_log10()

# effect on rich
ggplot(div_expl, aes(x = edge, y = Rich)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~taxa, scales = "free") 
  

ggplot(div_expl, aes(x = prox, y = Rich)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~taxa, scales = "free") 

ggplot(div_expl, aes(x = speccomp, y = Rich)) +
  geom_boxplot() +
  facet_wrap(~taxa, scales = "free") 

# effect on sha
ggplot(div_expl, aes(x = edge, y = Sha)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~taxa, scales = "free")

ggplot(div_expl, aes(x = prox, y = Sha)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~taxa, scales = "free")

ggplot(div_expl, aes(x = speccomp, y = Sha)) +
  geom_boxplot() +
  facet_wrap(~taxa, scales = "free")

# effect on Simp
ggplot(div_expl, aes(x = edge, y = Simp)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~taxa, scales = "free")

ggplot(div_expl, aes(x = prox, y = Simp)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~taxa, scales = "free")

ggplot(div_expl, aes(x = speccomp, y = Simp)) +
  geom_boxplot() +
  facet_wrap(~taxa, scales = "free")


m_a <- glm.nb(Abun ~ taxa * (edge + prox + speccomp),div_expl, contrasts = list(taxa="contr.sum"))
m_r <- lm(Rich ~ taxa * (edge + prox + speccomp),div_expl, contrasts = list(taxa="contr.sum"))
m_s <- lm(Sha ~ taxa * (edge + prox + speccomp),div_expl, contrasts = list(taxa="contr.sum"))
m_i <- lm(Simp ~ taxa * (edge + prox + speccomp),div_expl, contrasts = list(taxa="contr.sum"))

# brms model
m_div <- brm(Sha ~ edge + prox + speccomp + (edge + prox +speccomp | taxa), div_expl, family = "gaussian")

# some circlize stuff
fun_v <- read.csv("data/fun_v.csv",sep=" ")

## go through the indices
### abundance
cc_ab <- cor(cbind(func[,c(8,20,2,7,26,16,13,17)],div_all[,seq(5,33,4)]),use = "complete.obs")
cc_ab[upper.tri(cc_ab,diag=TRUE)] <- NA
cc_ab2 <- as.data.frame(cc_ab)
cc_ab2$From <- rownames(cc_ab2)
cc_ab3 <- gather(cc_ab2,"to","R",-From)
cc_ab3 <- subset(cc_ab3,!is.na(cc_ab3$R))

## richness
cc_ri <- cor(cbind(func[,c(8,20,2,7,26,16,13,17)],div_all[,seq(2,33,4)]),use = "complete.obs")
cc_ri[upper.tri(cc_ri,diag=TRUE)] <- NA
cc_ri2 <- as.data.frame(cc_ri)
cc_ri2$From <- rownames(cc_ri2)
cc_ri3 <- gather(cc_ri2,"to","R",-From)
cc_ri3 <- subset(cc_ri3,!is.na(cc_ri3$R))

## shannon
cc_sha <- cor(cbind(func[,c(8,20,2,7,26,16,13,17)],div_all[,seq(3,33,4)]),use = "complete.obs")
cc_sha[upper.tri(cc_sha,diag=TRUE)] <- NA
cc_sha2 <- as.data.frame(cc_sha)
cc_sha2$From <- rownames(cc_sha2)
cc_sha3 <- gather(cc_sha2,"to","R",-From)
cc_sha3 <- subset(cc_sha3,!is.na(cc_sha3$R))

## simpson
cc_sim <- cor(cbind(func[,c(8,20,2,7,26,16,13,17)],div_all[,seq(4,33,4)]),use = "complete.obs")
cc_sim[upper.tri(cc_sim,diag=TRUE)] <- NA
cc_sim2 <- as.data.frame(cc_sim)
cc_sim2$From <- rownames(cc_sim2)
cc_sim3 <- gather(cc_sim2,"to","R",-From)
cc_sim3 <- subset(cc_sim3,!is.na(cc_sim3$R))

# color indexing of the variables, change from ugly red/green to something nicer
grid.col_ab = c(Biomass = "red",Bird_smi = "red",C_stock="red",Decomp="red",frass="red",Herbivory="red",P_germ="red",Predation="red",Abun_veg="green",Abun_herb="green",Abun_cara="green",Abun_spider="green",Abun_iso="green",Abun_dip="green",Abun_bird="green",Abun_bat="green")
grid.col_ri = c(Biomass = "red",Bird_smi = "red",C_stock="red",Decomp="red",frass="red",Herbivory="red",P_germ="red",Predation="red",Rich_veg="green",Rich_herb="green",Rich_cara="green",Rich_spider="green",Rich_iso="green",Rich_dip="green",Rich_bird="green",Rich_bat="green")
grid.col_sim = c(Biomass = "red",Bird_smi = "red",C_stock="red",Decomp="red",frass="red",Herbivory="red",P_germ="red",Predation="red",Simp_veg="green",Simp_herb="green",Simp_cara="green",Simp_spider="green",Simp_iso="green",Simp_dip="green",Simp_bird="green",Simp_bat="green")
grid.col_sha = c(Biomass = "red",Bird_smi = "red",C_stock="red",Decomp="red",frass="red",Herbivory="red",P_germ="red",Predation="red",Sha_veg="green",Sha_herb="green",Sha_cara="green",Sha_spider="green",Sha_iso="green",Sha_dip="green",Sha_bird="green",Sha_bat="green")

# the plot
png("figures/correlation_all.png",width=1600,height=1200,pointsize = 18)
par(mfrow=c(2,2),cex=1.2)
circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,2),canvas.xlim=c(-1,1.5))
chordDiagram(cc_ab3,grid.col = grid.col_ab,col=colorRamp2(c(-0.33,0,0.33,0.67),c("#FDE725FF","white","#31688EFF","#440154FF")))
rect(1.1,0,1.2,0.8)
rasterImage(as.raster(matrix(c("#440154FF","#31688EFF","white","#FDE725FF"))),1.1,0,1.2,0.8)
text(x=1.1,y=1,labels="Correlation")
text(x=1.35,y=seq(0,0.8,length=4),labels = c(-0.33,0,0.33,0.66))
segments(rep(1.2,4),seq(0,0.8,length=4),rep(1.25,4),seq(0,0.8,length=4))
mtext("Correlation diagramm function - abundance")
circos.clear()

circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,2),canvas.xlim=c(-1,1.5))
chordDiagram(cc_ri3,grid.col = grid.col_ri,col=colorRamp2(c(-0.42,-0.06,0.3,0.66),c("#FDE725FF","white","#31688EFF","#440154FF")))
rect(1.1,0,1.2,0.8)
rasterImage(as.raster(matrix(c("#440154FF","#31688EFF","white","#FDE725FF"))),1.1,0,1.2,0.8)
text(x=1.1,y=1,labels="Correlation")
text(x=1.35,y=seq(0,0.8,length=4),labels = c(-0.42,-0.06,0.3,0.66))
segments(rep(1.2,4),seq(0,0.8,length=4),rep(1.25,4),seq(0,0.8,length=4))
mtext("Correlation diagramm function - richness")
circos.clear()

circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,2),canvas.xlim=c(-1,1.5))
chordDiagram(cc_sha3,grid.col = grid.col_sha,col=colorRamp2(c(-0.33,0,0.33,0.67),c("#FDE725FF","white","#31688EFF","#440154FF")))
rect(1.1,0,1.2,0.8)
rasterImage(as.raster(matrix(c("#440154FF","#31688EFF","white","#FDE725FF"))),1.1,0,1.2,0.8)
text(x=1.1,y=1,labels="Correlation")
text(x=1.35,y=seq(0,0.8,length=4),labels = c(-0.33,0,0.33,0.66))
segments(rep(1.2,4),seq(0,0.8,length=4),rep(1.25,4),seq(0,0.8,length=4))
mtext("Correlation diagramm function - shannon diversity")
circos.clear()

circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,2),canvas.xlim=c(-1,1.5))
chordDiagram(cc_sim3,grid.col = grid.col_sim,col=colorRamp2(c(-0.42,-0.06,0.3,0.67),c("#FDE725FF","white","#31688EFF","#440154FF")))
rect(1.1,0,1.2,0.8)
rasterImage(as.raster(matrix(c("#440154FF","#31688EFF","white","#FDE725FF"))),1.1,0,1.2,0.8)
text(x=1.1,y=1,labels="Correlation")
text(x=1.35,y=seq(0,0.8,length=4),labels = c(-0.42,-0.06,0.3,0.67))
segments(rep(1.2,4),seq(0,0.8,length=4),rep(1.25,4),seq(0,0.8,length=4))
mtext("Correlation diagramm function - simpson diversity")
circos.clear()

dev.off()


# some corr / cov analysis between residuals of the models
res_fun <- residuals(m_func,summary=FALSE)
res_div <- residuals(m_div,summary = FALSE)
res_abu <- residuals(m_abu,summary = FALSE)
# overall cov
hist(sapply(1:4000,function(i) cov(res_fun[i,],res_div[i,])))
hist(sapply(1:4000,function(i) cov(res_fun[i,],res_abu[i,])))
hist(sapply(1:4000,function(i) cov(res_abu[i,],res_div[i,])))

# variable - variable cov
aa_fun <- array(res_fun,dim=c(4000,53,8))
aa_div <- array(res_div,dim=c(4000,53,8))
aa_abu <- array(res_abu,dim=c(4000,53,8))
tt_fd <- sapply(1:4000,function(i) cor(cbind(aa_fun[i,,],aa_div[i,,])),simplify = FALSE)
tt_fa <- sapply(1:4000,function(i) cor(cbind(aa_fun[i,,],aa_abu[i,,])),simplify = FALSE)
tt_da <- sapply(1:4000,function(i) cor(cbind(aa_abu[i,,],aa_div[i,,])),simplify = FALSE)

# some formatting
## function - diversity
tmp <- array(as.numeric(unlist(tt_fd)),dim=c(16,16,4000),
             dimnames = list(c(unique(brm_fun$Function),paste0("Div_",unique(brm_div$group))),
                             c(unique(brm_fun$Function),paste0("Div_",unique(brm_div$group))),
                             paste0("Iter_",1:4000)))
res_cov <- adply(tmp,c(1,2),median) # might also get more than median

res_mat <- as.matrix(spread(res_cov,X2,V1)[,-1])
res_mat[upper.tri(res_mat,diag=TRUE)] <- NA
res_mat <- as.data.frame(res_mat)
res_mat$From <- c(unique(brm_dat_func$Function), paste0("Div_",unique(brm_div$group)))
res_fd <- gather(res_mat,"To","R",-From)
res_fd <- filter(res_fd,!is.na(R))
# set self cor to zero
res_fd$Type_from <- ifelse(res_fd$From %in% func_n, "function","diversity")
res_fd$Type_to <- ifelse(res_fd$To %in% func_n, "function","diversity")
res_fd$R <- ifelse(res_fd$Type_from == res_fd$Type_to,0,res_fd$R)
# get sign proba
si_fd <- adply(tmp,c(1,2),function(x) sum(x > 0) / 4000)
si_fd$V1 <- ifelse(si_fd$V1 < 0.5, 1 - si_fd$V1,si_fd$V1)
si_mat <- as.matrix(spread(si_fd,X2,V1)[,-1])
si_mat[upper.tri(si_mat,diag=TRUE)] <- NA
si_mat <- as.data.frame(si_mat)
si_mat$From <- c(unique(brm_fun$Function),paste0("Div_",unique(brm_div$group)))
sig_fd <- gather(si_mat,"To","P",-From)
sig_fd <- filter(sig_fd,!is.na(P))
# set self cor to zero
sig_fd$Type_from <- ifelse(sig_fd$From %in% func_n, "function","diversity")
sig_fd$Type_to <- ifelse(sig_fd$To %in% func_n, "function","diversity")
sig_fd$P <- ifelse(sig_fd$Type_from == sig_fd$Type_to,0,sig_fd$P)
# put back in original data frame
res_fd$P <- sig_fd$P
res_fd$P <- cut(res_fd$P,breaks = c(0,0.8,0.9,0.95,1),labels = c(0,3,2,1),include.lowest = TRUE) # turn into categoris
lty_fd <- as.numeric(levels(res_fd$P))[res_fd$P] # the lty

grid.col = c(Biomass = "red",Bird_smi = "red",C_stock="red",Decomp="red",frass="red",Herbivory="red",P_germ="red",Predation="red",
             Div_veg="green",Div_herb="green",Div_cara="green",Div_spider="green",Div_iso="green",Div_dip="green",Div_bird="green",Div_bat="green")

png("figures/rescor_fd.png",width=800,height=600,pointsize = 14)
circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,2),canvas.xlim=c(-1,1.5))
chordDiagram(res_fd[,1:3],grid.col = grid.col,col=colorRamp2(seq(-0.3,0.3,0.1),brewer.pal(7,"BrBG"))
             ,link.border = "black",link.lty = lty_fd)
rect(1.05,-0.45,1.66,1.05)
rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(7,"BrBG")[7:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=5),labels = c("-0.30","-0.15","0.00","0.15","0.30"))
segments(rep(1.2,4),seq(0.2,0.9,length=5),rep(1.25,4),seq(0.2,0.9,length=5))
mtext("Residual correlation diagramm function - diversity")
text(1.4,0,labels = "P(|Correlation| > 0)")
segments(rep(1.1,4),seq(-0.4,-0.1,length=4),rep(1.2,4),seq(-0.4,-0.1,length=4),lty = c(0,3,2,1))
text(x=rep(1.35,4),y=seq(-0.4,-0.1,length=4),labels = c("<0.80","0.80 - 0.90","0.90 - 0.95",">0.95"))
circos.clear()
dev.off()

## function - abundance
tmp <- array(as.numeric(unlist(tt_fa)),dim=c(16,16,4000),
             dimnames = list(c(unique(brm_fun$Function),paste0("Abu_",unique(brm_abu$group))),
                             c(unique(brm_fun$Function),paste0("Abu_",unique(brm_abu$group))),
                             paste0("Iter_",1:4000)))
res_cov <- adply(tmp,c(1,2),median) # might also get more than median

res_mat <- as.matrix(spread(res_cov,X2,V1)[,-1])
res_mat[upper.tri(res_mat,diag=TRUE)] <- NA
res_mat <- as.data.frame(res_mat)
res_mat$From <- c(unique(brm_dat_func$Function), paste0("Abu_",unique(brm_abu$group)))
res_fa <- gather(res_mat,"To","R",-From)
res_fa <- filter(res_fa,!is.na(R))
# set self cor to zero
res_fa$Type_from <- ifelse(res_fa$From %in% func_n, "function","abundance")
res_fa$Type_to <- ifelse(res_fa$To %in% func_n, "function","abundance")
res_fa$R <- ifelse(res_fa$Type_from == res_fa$Type_to,0,res_fa$R)
# get sign proba
si_fa <- adply(tmp,c(1,2),function(x) sum(x > 0) / 4000)
si_fa$V1 <- ifelse(si_fa$V1 < 0.5, 1 - si_fa$V1,si_fa$V1)
si_mat <- as.matrix(spread(si_fa,X2,V1)[,-1])
si_mat[upper.tri(si_mat,diag=TRUE)] <- NA
si_mat <- as.data.frame(si_mat)
si_mat$From <- c(unique(brm_fun$Function),paste0("Abu_",unique(brm_abu$group)))
sig_fa <- gather(si_mat,"To","P",-From)
sig_fa <- filter(sig_fa,!is.na(P))
# set self cor to zero
sig_fa$Type_from <- ifelse(sig_fa$From %in% func_n, "function","abundance")
sig_fa$Type_to <- ifelse(sig_fa$To %in% func_n, "function","abundance")
sig_fa$P <- ifelse(sig_fa$Type_from == sig_fa$Type_to,0,sig_fa$P)
# put back in original data frame
res_fa$P <- sig_fa$P
res_fa$P <- cut(res_fa$P,breaks = c(0,0.8,0.9,0.95,1),labels = c(0,3,2,1),include.lowest = TRUE) # turn into categoris
lty_fa <- as.numeric(levels(res_fa$P))[res_fa$P] # the lty

grid.col = c(Biomass = "red",Bird_smi = "red",C_stock="red",Decomp="red",frass="red",Herbivory="red",P_germ="red",Predation="red",
             Abu_veg="blue",Abu_herb="blue",Abu_cara="blue",Abu_spider="blue",Abu_iso="blue",Abu_dip="blue",Abu_bird="blue",Abu_bat="blue")

png("figures/rescor_fa.png",width=800,height=600,pointsize = 14)
circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,2),canvas.xlim=c(-1,1.5))
chordDiagram(res_fa[,1:3],grid.col = grid.col,col=colorRamp2(seq(-0.3,0.3,0.1),brewer.pal(7,"BrBG"))
             ,link.border = "black",link.lty = lty_fa)
rect(1.05,-0.45,1.66,1.05)
rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(7,"BrBG")[7:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=5),labels = c("-0.30","-0.15","0.00","0.15","0.30"))
segments(rep(1.2,4),seq(0.2,0.9,length=5),rep(1.25,4),seq(0.2,0.9,length=5))
mtext("Residual correlation diagramm function - abundance")
text(1.4,0,labels = "P(|Correlation| > 0)")
segments(rep(1.1,4),seq(-0.4,-0.1,length=4),rep(1.2,4),seq(-0.4,-0.1,length=4),lty = c(0,3,2,1))
text(x=rep(1.35,4),y=seq(-0.4,-0.1,length=4),labels = c("<0.80","0.80 - 0.90","0.90 - 0.95",">0.95"))
circos.clear()
dev.off()

## diversity - abundance
tmp <- array(as.numeric(unlist(tt_da)),dim=c(16,16,4000),
             dimnames = list(c(paste0("Abu_",unique(brm_abu$group)),paste0("Div_",unique(brm_div$group))),
                             c(paste0("Abu_",unique(brm_abu$group)),paste0("Div_",unique(brm_div$group))),
                             paste0("Iter_",1:4000)))
res_cov <- adply(tmp,c(1,2),median) # might also get more than median

res_mat <- as.matrix(spread(res_cov,X2,V1)[,-1])
res_mat[upper.tri(res_mat,diag=TRUE)] <- NA
res_mat <- as.data.frame(res_mat)
res_mat$From <- c(paste0("Abu_",unique(brm_abu$group)), paste0("Div_",unique(brm_div$group)))
res_da <- gather(res_mat,"To","R",-From)
res_da <- filter(res_da,!is.na(R))
# set self cor to zero
res_da$Type_from <- ifelse(grepl("Div_",res_da$From), "diversity","abundance")
res_da$Type_to <- ifelse(grepl("Div_",res_da$To), "diversity","abundance")
res_da$R <- ifelse(res_da$Type_from == res_da$Type_to,0,res_da$R)
# now get sign proba
si_da <- adply(tmp,c(1,2),function(x) sum(x > 0) / 4000)
si_da$V1 <- ifelse(si_da$V1 < 0.5, 1 - si_da$V1,si_da$V1)
si_mat <- as.matrix(spread(si_da,X2,V1)[,-1])
si_mat[upper.tri(si_mat,diag=TRUE)] <- NA
si_mat <- as.data.frame(si_mat)
si_mat$From <- c(paste0("Abu_",unique(brm_abu$group)), paste0("Div_",unique(brm_div$group)))
sig_da <- gather(si_mat,"To","P",-From)
sig_da <- filter(sig_da,!is.na(P))
# set self cor to zero
sig_da$Type_from <- ifelse(grepl("Div_",sig_da$From), "diversity","abundance")
sig_da$Type_to <- ifelse(grepl("Div_",sig_da$To), "diversity","abundance")
sig_da$P <- ifelse(sig_da$Type_from == sig_da$Type_to,0,sig_da$P)
# put back in original data frame
res_da$P <- sig_da$P
res_da$P <- cut(res_da$P,breaks = c(0,0.8,0.9,0.95,1),labels = c(0,3,2,1),include.lowest = TRUE) # turn into categoris
lty_da <- as.numeric(levels(res_da$P))[res_da$P] # the lty

grid.col = c(Abu_bat = "blue", Abu_bird = "blue", Abu_spider = "blue", Abu_herb = "blue", Abu_veg = "blue", Abu_iso = "blue", Abu_dip = "blue", Abu_cara = "blue",
             Div_bat = "green", Div_bird = "green", Div_spider = "green", Div_herb = "green", Div_veg = "green", Div_iso = "green", Div_dip = "green", Div_cara = "green")
  
png("figures/rescor_da.png",width=800,height=600,pointsize = 14)  
circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,2),canvas.xlim=c(-1,1.5))
chordDiagram(res_da[,1:3],grid.col = grid.col,col=colorRamp2(seq(-0.6,0.6,0.2),brewer.pal(7,"BrBG"))
             ,link.border = "black",link.lty = lty_da)
rect(1.05,-0.45,1.66,1.05)
rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(7,"BrBG")[7:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=5),labels = c("-0.60","-0.30","0.00","0.30","0.60"))
segments(rep(1.2,4),seq(0.2,0.9,length=5),rep(1.25,4),seq(0.2,0.9,length=5))
mtext("Residual correlation diagramm diversity - abundance")
text(1.4,0,labels = "P(|Correlation| > 0)")
segments(rep(1.1,4),seq(-0.4,-0.1,length=4),rep(1.2,4),seq(-0.4,-0.1,length=4),lty = c(0,3,2,1))
text(x=rep(1.35,4),y=seq(-0.4,-0.1,length=4),labels = c("<0.80","0.80 - 0.90","0.90 - 0.95",">0.95"))
circos.clear()
dev.off()

# all together !!!
res_all <- rbind(res_fd, res_fa, res_da)

grid.col = c(Abu_bat = "blue", Abu_bird = "blue", Abu_spider = "blue", Abu_herb = "blue", Abu_veg = "blue", Abu_iso = "blue", Abu_dip = "blue", Abu_cara = "blue",
             Div_bat = "green", Div_bird = "green", Div_spider = "green", Div_herb = "green", Div_veg = "green", Div_iso = "green", Div_dip = "green", Div_cara = "green",
             Biomass = "red", C_stock = "red", P_germ = "red", Decomp = "red", Herbivory = "red", Predation = "red", frass = "red", Bird_smi = "red")

lty_all <- as.numeric(levels(res_all$P))[res_all$P]

png("figures/rescor_all.png",width=800,height=600,pointsize = 14)
circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,rep(2,7),15,2),canvas.xlim = c(-1,1.5))
chordDiagram(res_all[,1:3], grid.col = grid.col, col=colorRamp2(seq(-0.6,0.6,0.2),brewer.pal(7,"BrBG")),
             link.lty = lty_all,link.border = "black")
rect(1.05,-0.45,1.66,1.05)
rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(7,"BrBG")[7:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=5),labels = c("-0.60","-0.30","0.00","0.30","0.60"))
segments(rep(1.2,4),seq(0.2,0.9,length=5),rep(1.25,4),seq(0.2,0.9,length=5))
mtext("Residual correlation diagramm function - diversity - abundance")
text(1.4,0,labels = "P(|Correlation| > 0)")
segments(rep(1.1,4),seq(-0.4,-0.1,length=4),rep(1.2,4),seq(-0.4,-0.1,length=4),lty = c(0,3,2,1))
text(x=rep(1.35,4),y=seq(-0.4,-0.1,length=4),labels = c("<0.80","0.80 - 0.90","0.90 - 0.95",">0.95"))
circos.clear()
dev.off()

# cool stuff now try to get sign prob for each corr coeff
tmp <- array(as.numeric(unlist(tt_da)),dim=c(16,16,4000),
             dimnames = list(c(paste0("Abu_",unique(brm_abu$group)),paste0("Div_",unique(brm_div$group))),
                             c(paste0("Abu_",unique(brm_abu$group)),paste0("Div_",unique(brm_div$group))),
                             paste0("Iter_",1:4000)))

########### scenario epxloration ########
# explore some scenario and its predicted impact on the different predictors
# for constant forest area, per plot

## 1. eradication of qrub (average changes from all -> fsyl_qrob, qrob_qrub -> qrob, fsyl_qrub -> fsyl)
newdata_fun <- expand.grid(speccomb = c("all","fsyl_qrob","qrob_qrub","qrob","fsyl_qrub","fsyl"), edge_std = 0, prox_std = 0, Function = unique(brm_fun$Function))
newdata_div <- expand.grid(speccomb = c("all","fsyl_qrob","qrob_qrub","qrob","fsyl_qrub","fsyl"), edge_std = 0, prox_std = 0, group = unique(brm_div$group))
newdata_abu <- expand.grid(speccomb = c("all","fsyl_qrob","qrob_qrub","qrob","fsyl_qrub","fsyl"), edge_std = 0, prox_std = 0, group = unique(brm_abu$group))

pred_fun <- posterior_linpred(m_fun, newdata = newdata_fun, summary = FALSE)
pred_div <- posterior_linpred(m_div, newdata = newdata_div, summary = FALSE)
pred_abu <- posterior_linpred(m_abu, newdata = newdata_abu, summary = FALSE)

### get percent change together with intervals for all functions
ch_fun <- adply(((pred_fun[,seq(2,48,2)] - pred_fun[,seq(1,48,2)]) / abs(pred_fun[,seq(1,48,2)])) * 100,2,quantile, probs = c(0.2,0.5,0.8))
ch_fun$Function <- rep(unique(brm_fun$Function),each=3)
ch_fun$Change <- rep(c("all -> fsyl_qrob","qrob_qrub -> qrob","fsyl_qrub -> fsyl"),8)
names(ch_fun)[2:4] <- c("LCI","Median","UCI")
ch_fun$X1 <- "Function"

ggplot(ch_fun,aes(x=Function,y=Median,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  facet_grid(~Change)

# for color coding
ch_fun2 <- adply(((pred_fun[,seq(2,48,2)] - pred_fun[,seq(1,48,2)]) / abs(pred_fun[,seq(1,48,2)])) * 100,2,function(x) c(sum(x > 0) / 4000,median(x)))
ch_fun2$V1 <- ifelse(ch_fun2$V1 < 0.5, 1 - ch_fun2$V1,ch_fun2$V1)
ch_fun2$group <- rep(unique(brm_fun$Function),each=3)
ch_fun2$Change <- rep(c("all -> fsyl_qrob","qrob_qrub -> qrob","fsyl_qrub -> fsyl"),8)

ch_div <- adply(((pred_div[,seq(2,48,2)] - pred_div[,seq(1,48,2)]) / abs(pred_div[,seq(1,48,2)])) * 100,2,quantile, probs = c(0.2,0.5,0.8))
ch_div$group <- rep(unique(brm_div$group),each=3)
ch_div$Change <- rep(c("all -> fsyl_qrob","qrob_qrub -> qrob","fsyl_qrub -> fsyl"),8)
names(ch_div)[2:4] <- c("LCI","Median","UCI")
ch_div$X1 <- "group"

ggplot(ch_div,aes(x=group,y=Median,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  facet_grid(~Change)
# for color coding
ch_div2 <- adply(((pred_div[,seq(2,48,2)] - pred_div[,seq(1,48,2)]) / abs(pred_div[,seq(1,48,2)])) * 100,2,function(x) c(sum(x > 0) / 4000,median(x)))
ch_div2$V1 <- ifelse(ch_div2$V1 < 0.5, 1 - ch_div2$V1,ch_div2$V1)
ch_div2$group <- rep(unique(brm_div$group),each=3)
ch_div2$Change <- rep(c("all -> fsyl_qrob","qrob_qrub -> qrob","fsyl_qrub -> fsyl"),8)

ch_abu <- adply(((pred_abu[,seq(2,48,2)] - pred_abu[,seq(1,48,2)]) / abs(pred_abu[,seq(1,48,2)])) * 100,2,quantile, probs = c(0.2,0.5,0.8))
ch_abu$group <- rep(unique(brm_abu$group),each=3)
ch_abu$Change <- rep(c("all -> fsyl_qrob","qrob_qrub -> qrob","fsyl_qrub -> fsyl"),8)
names(ch_abu)[2:4] <- c("LCI","Median","UCI")
ch_abu$X1 <- "group"

ggplot(ch_abu,aes(x=group,y=Median,ymin=LCI,ymax=UCI)) +
  geom_linerange() +
  geom_point() +
  facet_grid(~Change)
# for color coding
ch_abu2 <- adply(((pred_abu[,seq(2,48,2)] - pred_abu[,seq(1,48,2)]) / abs(pred_abu[,seq(1,48,2)])) * 100,2,function(x) c(sum(x > 0) / 4000,median(x)))
ch_abu2$V1 <- ifelse(ch_abu2$V1 < 0.5, 1 - ch_abu2$V1,ch_abu2$V1)
ch_abu2$group <- rep(unique(brm_abu$group),each=3)
ch_abu2$Change <- rep(c("all -> fsyl_qrob","qrob_qrub -> qrob","fsyl_qrub -> fsyl"),8)

# formatting for latex
ch_all2 <- rbind(ch_fun2[,3:5],ch_div2[,3:5],ch_abu2[,3:5])
ch_all2$V2 <- round(ch_all2$V2,2)
ch_all2$type <- rep(c("function","diversity","abundance"),each=8*3)

ch_dd <- dcast(ch_all2,group+type~Change,value.var = "V2")
ch_dd$type <- factor(ch_dd$type, levels = c("function","diversity","abundance"))
ch_dd <- arrange(ch_dd,type)


## 2. land sparing (constant composition, low edge, low proximity, compared to average, edge and average proximity)
newdata_fun <- expand.grid(speccomb = "all", edge_std = c(-2,0),prox_std = c(-2,0), Function = unique(brm_fun$Function))[c(seq(1,32,4),seq(4,32,4)),]
newdata_div <- expand.grid(speccomb = "all", edge_std = c(-2,0),prox_std = c(-2,0), group = unique(brm_div$group))[c(seq(1,32,4),seq(4,32,4)),]
newdata_abu <- expand.grid(speccomb = "all", edge_std = c(-2,0),prox_std = c(-2,0), group = unique(brm_div$group))[c(seq(1,32,4),seq(4,32,4)),]

pred_fun <- posterior_linpred(m_fun, newdata = newdata_fun, summary = FALSE)
pred_div <- posterior_linpred(m_div, newdata = newdata_div, summary = FALSE)
pred_abu <- posterior_linpred(m_abu, newdata = newdata_abu, summary = FALSE)

## for the function turn back in the original scale
mm_fun <- mean

## percent change this time directly with proba
ch_fun <- adply(((pred_fun[,1:8] - pred_fun[,9:16]) / abs(pred_fun[,9:16])) * 100,2,function(x) c(sum(x > 0) / 4000, median(x)))

## 3. land sharing (constant composition, high edge, average proximity, compared to average edge)


########## observed (raw) correlations between variables ########
dat_fd <- cbind(func[,c("C_stock","Biomass","P_germ","Herbivory","Decomp","Predation","frass","Bird_smi")],
                div_all[,c("Sha_veg","Sha_herb","Sha_iso","Sha_dip","Sha_cara","Sha_spider","Sha_bird","Sha_bat")])

cc_fd <- cor(dat_fd,use = "complete.obs")
cc_fd[upper.tri(cc_fd,diag=TRUE)] <- NA
cc_fd <- as.data.frame(cc_fd)
cc_fd$From <- rownames(cc_fd)
cc_fd %>%
  gather("To","R",-From) %>%
  filter(!is.na(R)) %>%
  mutate(Type_from = ifelse(From %in% unique(brm_fun$Function), "function","diversity")) %>%
  mutate(Type_to = ifelse(To %in% unique(brm_fun$Function), "function","diversity")) %>%
  mutate(R = ifelse(Type_from == Type_to,0,R)) -> cc_fd2


grid.col = c(Biomass = "red",Bird_smi = "red",C_stock="red",Decomp="red",frass="red",Herbivory="red",P_germ="red",Predation="red",
             Sha_veg="blue",Sha_herb="blue",Sha_cara="blue",Sha_spider="blue",Sha_iso="blue",Sha_dip="blue",Sha_bird="blue",Sha_bat="blue")

png("figures/observed_fd.png",width=1000,height=800,pointsize = 18)
circos.par(gap.after=c(rep(2,6),15,rep(2,7),15,2),canvas.xlim=c(-1,1.5))
chordDiagram(cc_fd2[,1:3],col=colorRamp2(seq(-0.35,0.35,0.1),brewer.pal(8,"BrBG")),grid.col = grid.col)
rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(8,"BrBG")[8:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=4),labels = c("-0.35","-0.15","0.15","0.35"))
segments(rep(1.2,4),seq(0.2,0.9,length=4),rep(1.25,4),seq(0.2,0.9,length=4))
mtext("Observed correlation diagramm function - diversity")
circos.clear()
dev.off()

# now do the same for function only
cc_f <- cor(func[,c("C_stock","Biomass","P_germ","Herbivory","Decomp","Predation","frass","Bird_smi")],use = "complete.obs")
cc_f[upper.tri(cc_f,diag=TRUE)] <- NA
cc_f <- as.data.frame(cc_f)
cc_f$From <- rownames(cc_f)
cc_f %>%
  gather("To","R",-From) %>%
  filter(!is.na(R)) -> cc_f2

grid.col = c(Biomass = "red",Bird_smi = "red",C_stock="red",Decomp="red",frass="red",Herbivory="red",P_germ="red",Predation="red")

png("figures/observed_f.png",width=1000,height=800,pointsize = 22)
circos.par(canvas.xlim=c(-1,1.5))
chordDiagram(cc_f2[,1:3],col=colorRamp2(seq(-0.67,0.67,length=5),brewer.pal(5,"BrBG")),grid.col = grid.col)
rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(5,"BrBG")[5:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=5),labels = c("-0.67","-0.33","0.00","0.33","0.67"))
segments(rep(1.2,4),seq(0.2,0.9,length=5),rep(1.25,4),seq(0.2,0.9,length=5))
mtext("Observed correlation diagramm function - function")
circos.clear()
dev.off()

# now do the same for diversity
cc_d <- cor(div_all[,c("Sha_veg","Sha_herb","Sha_iso","Sha_dip","Sha_cara","Sha_spider","Sha_bird","Sha_bat")],use = "complete.obs")
cc_d[upper.tri(cc_d,diag=TRUE)] <- NA
cc_d <- as.data.frame(cc_d)
cc_d$From <- rownames(cc_d)
cc_d %>%
  gather("To","R",-From) %>%
  filter(!is.na(R)) -> cc_d2

grid.col = c( Sha_veg="blue",Sha_herb="blue",Sha_cara="blue",Sha_spider="blue",Sha_iso="blue",Sha_dip="blue",Sha_bird="blue",Sha_bat="blue")

png("figures/observed_d.png",width=1000,height=800,pointsize = 22)
circos.par(canvas.xlim=c(-1,1.5))
chordDiagram(cc_d2[,1:3],col=colorRamp2(seq(-0.35,0.35,length=5),brewer.pal(5,"BrBG")),grid.col = grid.col)
rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(5,"BrBG")[5:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=5),labels = c("-0.35","-0.17","0.00","0.17","0.35"))
segments(rep(1.2,4),seq(0.2,0.9,length=5),rep(1.25,4),seq(0.2,0.9,length=5))
mtext("Observed correlation diagramm diversity - diversity")
circos.clear()
dev.off()

####### start brms multivariate analysis #####
brm_mult <- cbind(func[,c(1,2,7,8,13,16,17,20,26)],div_all[,seq(3,33,4)],pred_dat[,-1])
brm_mult <- mutate_at(brm_mult,2:19,scale)
brm_mult[is.na(brm_mult)] <- 0

m_mult <- brm(formula = mvbind(C_stock, Decomp, Biomass, P_germ, Herbivory, Predation, Bird_smi, frass, Sha_veg, Sha_cara, Sha_spider, Sha_iso, Sha_dip, Sha_bird, Sha_bat) ~
                frag_prox + plot_edge100m + speccomb, data = brm_mult)

m_int <- brm(formula = mvbind(C_stock, Decomp, Biomass, P_germ, Herbivory, Predation, Bird_smi, frass, Sha_veg, Sha_cara, Sha_spider, Sha_iso, Sha_dip, Sha_bird, Sha_bat) ~
               (frag_prox + plot_edge100m) * speccomb, data = brm_mult)

m_null <- brm(formula = mvbind(C_stock, Decomp, Biomass, P_germ, Herbivory, Predation, Bird_smi, frass, Sha_veg, Sha_cara, Sha_spider, Sha_iso, Sha_dip, Sha_bird, Sha_bat) ~
               1, data = brm_mult)
p_null <- posterior_samples(m_null)
p_null[,31:135] %>%
  gather("corr","value") %>%
  group_by(corr) %>%
  summarise(R = median(value), P = sum(value > 0) / 4000) %>%
  mutate(P = ifelse(P < 0.5, 1 - P, P)) %>%
  arrange(desc(P)) %>%
  separate(corr,c("drop","From","To")) %>%
  select(-drop) %>%
  mutate(Type_from = ifelse(From %in% unique(brm_fun$Function),"function","diversity")) %>%
  mutate(Type_to = ifelse(To %in% unique(brm_fun$Function),"function","diversity")) -> null_dd


p_mult <- posterior_samples(m_mult)
p_mult[,151:255] %>%
  gather("corr","value") %>%
  group_by(corr) %>%
  summarise(R_mult = median(value), P_mult = sum(value > 0) / 4000) %>%
  mutate(P_mult = ifelse(P_mult < 0.5, 1 - P_mult, P_mult)) %>%
  arrange(desc(P_mult)) %>%
  separate(corr,c("drop","From","To")) %>%
  select(-drop) %>%
  mutate(Type_from = ifelse(From %in% unique(brm_fun$Function),"function","diversity")) %>%
  mutate(Type_to = ifelse(To %in% unique(brm_fun$Function),"function","diversity")) -> mult_dd

dd_all <- left_join(null_dd,mult_dd,by=c("From","To"))
dd_all$change <- "none"
dd_all$change <- ifelse(dd_all$P > 0.8 & dd_all$P_mult < 0.8, "extrinsic",dd_all$change)
dd_all$change <- ifelse(dd_all$P > 0.8 & dd_all$P_mult > 0.8, "intrinsic",dd_all$change)
dd_all$change <- ifelse(dd_all$P < 0.8 & dd_all$P_mult > 0.8, "appearing",dd_all$change)

dd_all$change_type <- paste(dd_all$Type_from.x, dd_all$Type_to.x,sep="_")
dd_all$change_type[dd_all$change_type == "diversity_function"] <- "function_diversity"
dd_all %>%
  group_by(change_type) %>%
  summarise(Ex_n = sum(change == "extrinsic"),
            In_n = sum(change == "intrinsic"),
            Ap_n = sum(change == "appearing")) %>%
  gather("direction","number",-change_type) %>%
  group_by(change_type) %>%
  mutate(prop = number / sum(number)) %>%
  arrange(change_type) -> dd_2

gg_cor <- ggplot(dd_2,aes(x=direction,y=prop,fill=change_type)) +
  geom_bar(stat="identity",position="dodge") +
  scale_x_discrete(labels = c("Appearing correlation","Extrinsic correlation","Intrinsic correlation"), name = "Type of correlation") +
  labs(y = "Proportion of correlation changes")

ggsave("figures/correlation_shift.png",gg_cor)


############# multivariate all together in the joy ##########

# a helper function
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

names(div_all)[2:33] <- gsub("_",".",names(div_all)[2:33])
# create the dataset
div_all %>%
  select(id_plot, contains("Sha")) %>%
  mutate_at(2:9, scale_this) %>%
  arrange(id_plot) -> div_dd

func %>%
  select(id_plot,C_stock,Decomp,Biomass,P_germ,Herbivory,Predation,Bird_smi,frass) %>%
  rename(P.germ = P_germ, C.stock = C_stock, Bird.smi = Bird_smi) %>%
  mutate_at(2:9, scale_this) %>%
  left_join(div_dd, by = "id_plot") %>%
  left_join(pred_dat, by = "id_plot") %>%
  mutate(prox.std = scale_this(frag_prox), edge.std = scale_this(plot_edge100m), dens.std = scale_this(total_all_ba)) %>%
  mutate(speccomb = gsub("_",".",speccomb))-> brm_dat
  
brm_dat[is.na(brm_dat)] <- 0 # give average values to NA rows


# fit the model
m_mult <- brm(mvbind(C.stock,Decomp,Biomass,P.germ,Herbivory,Predation,Bird.smi,frass,Sha.veg,Sha.herb,Sha.cara,Sha.spider,Sha.iso,Sha.dip,Sha.bird,Sha.bat)
              ~ speccomb + prox.std + edge.std + dens.std,
              prior = set_prior("normal(0, 2)", class = "b"), data = brm_dat)
# elpd_kfold  -1223.9 (SE: 32.7)

# interaction try out
m_int <- brm(mvbind(C.stock,Decomp,Biomass,P.germ,Herbivory,Predation,Bird.smi,frass,Sha.veg,Sha.herb,Sha.cara,Sha.spider,Sha.iso,Sha.dip,Sha.bird,Sha.bat)
              ~ speccomb * (prox.std + edge.std),
              prior = set_prior("normal(0, 2)", class = "b"), data = brm_dat)

# or load it from a file
m_mult <- readRDS("model/m_mult.rds")


# 0. posterior predictive checks
nn <- names(brm_dat)[2:17]
nn <- gsub("\\.","",nn)
gg <- list()

for(i in 1:16){
 gg[[i]] <- pp_check(m_mult, resp = nn[i], nsamples = 100) + labs(title = nn[i])
}

ag <- arrangeGrob(grobs = gg)

# 1. look at the effects
cc <- as.data.frame(fixef(m_mult))
## 1.1 fragmentation effects (prox + edge)
cc$id <- rownames(cc)

cc %>%
  separate(id,c("variable","parameter"),sep = "_") %>%
  arrange(variable) %>%
  mutate(type = rep(c("function","diversity"),each = 9 * 8)) %>%
  filter(parameter %in% c("edge.std","prox.std")) -> frag_dd

name_nice <- c("Tree biomass", "Bird biomass", "Carbon stock", "Decomposition", "Insect biomass", "Herbivory", "Seedling success", "Predation", "Bat", "Bird", "Carabid", "Diplopod", "Herbivore", "Isopod", "Spider", "Vegetation")

frag_dd$variable <- rep(name_nice, each = 2)
frag_dd$parameter <- rep(c("proximity","edge"),times = 16)
frag_dd$X <- rep(1:16,each=2) + rep(c(-0.1,0.1),times=16)


### get average response
frag_dd %>%
  group_by(parameter, type) %>%
  summarise(Estimate = mean(Estimate)) -> frag_avg

gg_frag <- ggplot(frag_dd, aes(x=X,y=Estimate, ymin = Q2.5, ymax = Q97.5,color = parameter)) +
  geom_linerange() +
  geom_point() + 
  coord_flip() +
  facet_wrap(~type, scales = "free_y") +
  labs(x="", y = "Parameter value (95% Credible interval)") +
  scale_x_continuous(breaks = 1:16, labels = name_nice) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(text = element_text(size = 20))


ggsave("figures/frag_eff.png",gg_frag, width = 12, height = 8)

## 1.2 composition effect
newdat <- expand.grid(speccomb = c("all","fsyl","fsyl.qrob","qrob.qrub","qrob","fsyl.qrub","qrub"),prox.std=0,edge.std=0)

comp_eff <- fitted(m_mult, newdata = newdat)
comp_eff <- adply(comp_eff,3,identity)
comp_eff$speccomb <- rep(newdat$speccomb, times = 16)
comp_eff$type <- rep(c("function","diversity"),each = 8*7)
comp_eff$X1 <- as.character(comp_eff$X1)
comp_eff <- arrange(comp_eff, X1)
comp_eff$X <- rep(1:16, each = 7) + rep(c(0.3,-0.3,0,0.2,-0.2,0.1,-0.1),times = 16)
comp_eff$speccomb <- factor(comp_eff$speccomb, levels = c("fsyl","qrob","qrub","fsyl.qrob","fsyl.qrub","qrob.qrub","all"))

ll <- data.frame(type = rep(c("function","diversity"),each=8), line = seq(1.5,16.5,by = 1))

gg_comp <- ggplot(comp_eff,aes(x=X,y=Estimate, ymin = Q2.5, ymax = Q97.5, color = speccomb)) +
  geom_linerange() +
  geom_point() +
  facet_wrap(~type, scales = "free_y",ncol = 2) +
  coord_flip() +
  scale_x_continuous(breaks = 1:16, labels = name_nice, name = "") +
  geom_vline(data =ll, aes(xintercept = line), linetype = "dashed") +
  labs(y = "Parameter value (95% Credible interval)") +
  scale_color_discrete(name = "Tree species\ncomposition") +
  theme(text = element_text(size = 20))


ggsave("figures/comp_eff.png", width = 12, height = 8)

## 1.3 R-square
r2 <- as.data.frame(bayes_R2(m_mult))
r2$variable <- substr(rownames(r2),3,nchar(rownames(r2)))
r2$type <- rep(c("function","diversity"),each=8)
r2 %>%
  arrange(variable) %>%
  mutate(variable = name_nice) -> r2

gg_r2 <- ggplot(r2,aes(x=variable,y=Estimate,ymin=Q2.5,ymax=Q97.5)) +
  geom_linerange() +
  geom_point() +
  facet_wrap(~type,scales = "free_y") +
  coord_flip() +
  labs(x="", y = "Parameter value (95% Credible interval)")

ggsave("figures/r2_fig.png", gg_r2, width = 12, height = 8)

# 2. make the scenarios

## 2.1 eradication of qrub
### the newdat
newdat <- data.frame(speccomb = c("all","fsyl.qrob","fsyl.qrub","fsyl","qrob.qrub","qrob"), prox.std = 0, edge.std = 0)

pred <- posterior_linpred(m_mult, newdata = newdat)
### put back in the original scale
func %>%
  select(id_plot,C_stock,Decomp,Biomass,P_germ,Herbivory,Predation,Bird_smi,frass) %>%
  rename(P.germ = P_germ, C.stock = C_stock, Bird.smi = Bird_smi) %>%
  gather("variable","value",-id_plot) %>%
  group_by(variable) %>%
  summarise(mean = mean(value,na.rm=TRUE),sd = sd(value,na.rm=TRUE)) -> mean_fun
  
div_all %>%
  select(id_plot, contains("Sha")) %>%
  gather("variable","value",-id_plot) %>%
  group_by(variable) %>%
  summarise(mean = mean(value,na.rm=TRUE),sd = sd(value,na.rm=TRUE)) -> mean_div

mean_all <- rbind(mean_fun, mean_div)
mean_all$variable <- gsub("\\.","",mean_all$variable)

for(i in 1:16){
  dd <- dimnames(pred)[[3]][i]
  m <- as.numeric(mean_all[mean_all$variable == dd, "mean"])
  sd <- as.numeric(mean_all[mean_all$variable == dd, "sd"])
  pred[,,i] <- pred[,,i] * sd + m  
}
ch_out <- NULL
for(i in 1:16){
  dd <- dimnames(pred)[[3]][i]
  ch_tmp <- ((pred[,seq(2,6,2),i] - pred[,seq(1,6,2),i]) / abs(pred[,seq(1,6,2),i])) * 100
  ch_tmp <- adply(ch_tmp, 2, function(x) c(sum(x > 0) / 4000, median(x))) 
  ch_tmp$change <- c("all -> fsyl_qrob","fsyl_qrub -> fsyl", "qrob_qrub -> qrob")
  ch_tmp$variable <- dd
  names(ch_tmp)[2:3] <- c("prob","percent")
  ch_tmp$prob <- ifelse(ch_tmp$prob < 0.5, 1 - ch_tmp$prob, ch_tmp$prob)
  ch_out <- rbind(ch_out, ch_tmp)
}

# some formatting
ch_out$percent <- round(ch_out$percent, 2)
ch_dd <- dcast(ch_out, variable ~ change, value.var = "percent")
xtable(ch_dd)


## 2.2 land sparing, low edge, low proximity
newdat <- expand.grid(speccomb = "all", edge.std = c(-2,0),prox.std = c(-2,0))[c(1,4),]
### prediction
pred <- posterior_linpred(m_mult, newdata = newdat)

for(i in 1:16){
  dd <- dimnames(pred)[[3]][i]
  m <- as.numeric(mean_all[mean_all$variable == dd, "mean"])
  sd <- as.numeric(mean_all[mean_all$variable == dd, "sd"])
  pred[,,i] <- pred[,,i] * sd + m  
}
ch_out <- NULL
for(i in 1:16){
  dd <- dimnames(pred)[[3]][i]
  ch_tmp <- ((pred[,1,i] - pred[,2,i]) / abs(pred[,2,i])) * 100
  ch_tmp <- data.frame(prob = sum(ch_tmp > 0) / 4000, percent = median(ch_tmp))
  ch_tmp$variable <- dd
  ch_tmp$prob <- ifelse(ch_tmp$prob < 0.5, 1 - ch_tmp$prob, ch_tmp$prob)
  ch_out <- rbind(ch_out, ch_tmp)
}

# some formatting
ch_out$percent <- round(ch_out$percent, 2)

xtable(ch_dd)



## 2.3 land sharing, high edge, average proximity
newdat <- expand.grid(speccomb = "all", edge.std = c(2,0),prox.std = 0)
### prediction
pred <- posterior_linpred(m_mult, newdata = newdat)

for(i in 1:16){
  dd <- dimnames(pred)[[3]][i]
  m <- as.numeric(mean_all[mean_all$variable == dd, "mean"])
  sd <- as.numeric(mean_all[mean_all$variable == dd, "sd"])
  pred[,,i] <- pred[,,i] * sd + m  
}
ch_sha <- NULL
for(i in 1:16){
  dd <- dimnames(pred)[[3]][i]
  ch_tmp <- ((pred[,1,i] - pred[,2,i]) / abs(pred[,2,i])) * 100
  ch_tmp <- data.frame(prob = sum(ch_tmp > 0) / 4000, percent = median(ch_tmp))
  ch_tmp$variable <- dd
  ch_tmp$prob <- ifelse(ch_tmp$prob < 0.5, 1 - ch_tmp$prob, ch_tmp$prob)
  ch_sha <- rbind(ch_sha, ch_tmp)
}

# some formatting
ch_sha$percent <- round(ch_sha$percent, 2)


ch_land <- data.frame(variable = ch_sha$variable, sparing = ch_out$percent, sharing = ch_sha$percent)
ch_land <- arrange(ch_land, variable)
xtable(ch_land)




# 3. look at the trade-offs

## 3.1 fit a null model
m_null <- brm(formula = mvbind(C.stock, Decomp, Biomass, P.germ, Herbivory, Predation, Bird.smi, frass, Sha.veg, Sha.herb, Sha.cara, Sha.spider, Sha.iso, Sha.dip, Sha.bird, Sha.bat) ~
                1, data = brm_dat)


p_null <- posterior_samples(m_null)
p_null[,33:152] %>%
  gather("corr","value") %>%
  group_by(corr) %>%
  summarise(R = median(value), P = sum(value > 0) / 4000) %>%
  mutate(P = ifelse(P < 0.5, 1 - P, P)) %>%
  arrange(desc(P)) %>%
  separate(corr,c("drop","From","To")) %>%
  select(-drop) %>%
  mutate(Type_from = ifelse(From %in% c("Cstock", "Decomp", "Biomass", "Pgerm", "Herbivory", "Predation", "Birdsmi", "frass"),"function","diversity")) %>%
  mutate(Type_to = ifelse(To %in% c("Cstock", "Decomp", "Biomass", "Pgerm", "Herbivory", "Predation", "Birdsmi", "frass"),"function","diversity")) -> null_dd

null_dd %>%
  filter(Type_from != Type_to) %>%
  mutate(From = factor(From), To = factor(To)) -> null_fd

null_fd$From <- revalue(null_fd$From, c("Biomass"="TBiom","Birdsmi"="BBiom","frass"="Insect biomass"))
null_fd$To <- revalue(null_fd$To, c("Shaveg"="Vegetation","Shaherb"="Herbivore","Shaiso"="Isopod","Shadip"="Dipl.","Shaspider"="Spider","Shacara"="Carabid","Shabird"="Bird","Shabat"="Bat"))
# reproduce figure from Felipe-Lucia on average levels of trade-offs
## the observed correlation
obs_r <- cor(brm_dat[,2:17])
obs_r[lower.tri(obs_r, diag=TRUE)] <- NA
obs_r <- melt(obs_r)
obs_r <- subset(obs_r, !is.na(value))
obs_r$col <- "a"
obs_r$col[obs_r$value > 0] <- "b"
obs_r$lab <- "Observed"

obs_r %>%
  group_by(col) %>%
  summarise(value = mean(value)) -> obs_avg


# grab the rescor post median
pp <- posterior_samples(m_mult)

rescor <- pp[grep("rescor",names(pp))]
rescor %>%
  gather(type, value) %>%
  group_by(type) %>%
  summarise(R = median(value), P = sum(value > 0) / n()) %>%
  mutate(P = ifelse(P < 0.5, 1 - P, P)) %>%
  mutate(Intrisic = ifelse(P > 0.9, "Y","N")) -> m_dd

# the observed corr
pp <- posterior_samples(m_null)

rescor <- pp[grep("rescor",names(pp))]
rescor %>%
  gather(type, value) %>%
  group_by(type) %>%
  summarise(R = median(value), P = sum(value > 0) / n()) %>%
  mutate(P = ifelse(P < 0.5, 1 - P, P)) %>%
  mutate(Present = ifelse(P > 0.9, "Y","N")) %>%
  left_join(m_dd[,c(1,4)], by = "type") %>%
  mutate(Res = paste(Present, Intrisic, sep = "_")) %>%
  separate(type, c(NA,"From","To"), sep  = "__") %>%
  mutate(From = factor(From, levels = ord), To = factor(To, levels = ord)) %>%
  mutate(label = ifelse(Res == "Y_Y",sprintf("underline(%s)",round(R,2)),round(R,2))) %>%
  mutate(label = ifelse(Res == "Y_N", sprintf("italic('%s')",round(R,2)), label)) %>%
  mutate(size = ifelse(Res %in% c("Y_N","Y_Y"), 2, 1)) %>%
  mutate(type_f = ifelse(From %in% fun, "function","diversity"),
         type_t = ifelse(To %in% fun, "function","diversity")) %>%
  unite("type",c(type_f,type_t),sep="_") %>%
  mutate(color = colorRamp2(c(-0.55,0,0.55),c("red","white","blue"))(R)) %>%
  mutate(From = revalue(From, c(Cstock = "C stock", Decomp = "Decomposition", Biomass = "Tree biomass", Pgerm = "Regeneration", Birdsmi = "Bird biomass", frass = "Insect biomass",
                                Shaveg = "Vegetation", Shaherb = "Herbivore", Shacara = "Carabid", Shaspider = "Spider", Shaiso = "Isopod", Shadip = "Diplopod", Shabird = "Bird", Shabat = "Bat"))) %>%
  mutate(To = revalue(To, c(Cstock = "C stock", Decomp = "Decomposition", Biomass = "Tree biomass", Pgerm = "Regeneration", Birdsmi = "Bird biomass", frass = "Insect biomass",
                                Shaveg = "Vegetation", Shaherb = "Herbivore", Shacara = "Carabid", Shaspider = "Spider", Shaiso = "Isopod", Shadip = "Diplopod", Shabird = "Bird", Shabat = "Bat"))) %>%
  mutate(type = revalue(type, c(diversity_diversity = "a) Diversity - Diversity", function_function = "b) Function - Function", function_diversity = "c) Function - Diversity"))) -> obs_dd


ord <- c("Cstock","Decomp", "Biomass","Pgerm","Herbivory","Predation","Birdsmi","frass",
         "Shaveg","Shaherb","Shacara","Shaspider","Shaiso","Shadip","Shabird","Shabat")

fun <- c("Cstock","Decomp", "Biomass","Pgerm","Herbivory","Predation","Birdsmi","frass")
div <- c("Shaveg","Shaherb","Shacara","Shaspider","Shaiso","Shadip","Shabird","Shabat")

obs_dd %>%
  

obs_dd$color <- colorRamp2(c(-0.55,0,0.55),c("red","white","blue"))(obs_dd$R)

gg_cor <- ggplot(obs_dd,aes(x=From,y=To,fill=color)) +
  geom_tile() +
  geom_text(aes(label = label, size = size), parse = TRUE, show.legend = FALSE) +
  scale_fill_identity() +
  facet_wrap(~type,scales = "free", ncol = 1) +
  scale_size(range = c(2,4)) +
  labs(x = "", y = "") +
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), axis.ticks = element_blank()) 


ggsave("figures/tradeoffs.png",gg_cor,width=7, height = 10)


mutate(col = ifelse(value < 0, "a", "b"), lab = "Composition\nfragm. removed") -> m_dd

  

# put together
r_all <- rbind(m_dd[,2:4], obs_r[,3:5])

ggplot(r_all, aes(x=lab, y= value, color = col)) +
  geom_point() +
  geom_segment(aes(x = 1.8, xend = 2.2, y = -0.116, yend = -0.116), color = "red") +
  geom_segment(aes(x = 1.8, xend = 2.2, y = 0.145, yend = 0.145), color = "blue")




grid.col = c(rep("salmon",8),rep("palegreen",8))
lty_fd <- cut(null_fd$P, breaks = c(0,0.8,0.9,0.95,1),labels = c(0,3,2,1), include.lowest = TRUE)
lty_fd <-  as.numeric(levels(lty_fd))[lty_fd]
lty_lwd <- ifelse(lty_fd == 1, 3, ifelse(lty_fd == 2, 1, ifelse(lty_fd == 3, 0.5, 0.1)))


png("figures/observed_fd.png",width=1000,height=800,pointsize = 18)
circos.par(gap.after=rep(2,16),canvas.xlim=c(-1,1.6))

chordDiagram(null_fd[,1:3], grid.col = grid.col, 
             annotationTrack = c("grid", "axis"), annotationTrackHeight = uh(9, "mm"),
             col=colorRamp2(seq(-0.26,0.26,length=8),brewer.pal(8,"BrBG")),
             link.lty = lty_fd, link.border = "black", link.lwd = lty_lwd)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black",cex=0.8)
}

rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(8,"BrBG")[8:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=4),labels = c("-0.26","-0.09","0.09","0.26"))
segments(rep(1.2,4),seq(0.2,0.9,length=4),rep(1.25,4),seq(0.2,0.9,length=4))
mtext("Observed correlation diagramm function - diversity")
text(1.4,0,labels = "P(|Correlation| > 0)")
segments(rep(1.1,4),seq(-0.4,-0.1,length=4),rep(1.2,4),seq(-0.4,-0.1,length=4),lty = c(0,3,2,1),lwd = c(0.1, 0.5, 1, 3))
text(x=rep(1.35,4),y=seq(-0.4,-0.1,length=4),labels = c("<0.80","0.80 - 0.90","0.90 - 0.95",">0.95"))
circos.clear()
dev.off()

# now for function - function
null_dd %>%
  filter(Type_from == "function" & Type_to == "function") %>%
  mutate(From = factor(From), To = factor(To)) -> null_f

null_f$From <- revalue(null_f$From, c("Biomass"="Tree biomass", "Decomp"="Decomposition","Cstock"="Carbon stock","Pgerm"="Seedling success","Birdsmi"="Bird biomass","frass"="Insect biomass"))
null_f$To <- revalue(null_f$To, c("Biomass"="Tree biomass", "Decomp"="Decomposition","Cstock"="Carbon stock","Pgerm"="Seedling success","Birdsmi"="Bird biomass","frass"="Insect biomass"))

grid.col = c(rep("salmon",8))
lty_fd <- cut(null_f$P, breaks = c(0,0.8,0.9,0.95,1),labels = c(0,3,2,1), include.lowest = TRUE)
lty_fd <-  as.numeric(levels(lty_fd))[lty_fd]
lty_lwd <- ifelse(lty_fd == 1, 3, ifelse(lty_fd == 2, 1, ifelse(lty_fd == 3, 0.5, 0.1)))


png("figures/observed_f.png",width=1000,height=800,pointsize = 18)
par(cex = 1.3)
circos.par(canvas.xlim=c(-1,1.6))
chordDiagram(null_f[,1:3], grid.col = grid.col, 
             annotationTrack = c("grid", "axis"), annotationTrackHeight = uh(9, "mm"),
             col=colorRamp2(seq(-0.54,0.54,length=8),brewer.pal(8,"BrBG")),
             link.lty = lty_fd, link.border = "black", link.lwd = lty_lwd)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black")
}

rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(8,"BrBG")[8:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=4),labels = c("-0.54","-0.18","0.18","0.54"))
segments(rep(1.2,4),seq(0.2,0.9,length=4),rep(1.25,4),seq(0.2,0.9,length=4))
mtext("Observed correlation diagramm function - function")
text(1.4,0,labels = "P(|Correlation| > 0)")
segments(rep(1.1,4),seq(-0.4,-0.1,length=4),rep(1.2,4),seq(-0.4,-0.1,length=4),lty = c(0,3,2,1),lwd = c(0.1, 0.5, 1, 3))
text(x=rep(1.35,4),y=seq(-0.4,-0.1,length=4),labels = c("<0.80","0.80 - 0.90","0.90 - 0.95",">0.95"))
circos.clear()
dev.off()

# now for diversity - diversity
null_dd %>%
  filter(Type_from == "diversity" & Type_to == "diversity") %>%
  mutate(From = factor(From), To = factor(To)) -> null_d

null_d$To <- revalue(null_d$To, c("Shaveg"="Vegetation","Shaherb"="Herbivore","Shaiso"="Isopod","Shadip"="Diplopod","Shaspider"="Spider","Shacara"="Carabid","Shabird"="Bird","Shabat"="Bat"))
null_d$From <- revalue(null_d$From, c("Shaveg"="Vegetation","Shaherb"="Herbivore","Shaiso"="Isopod","Shadip"="Diplopod","Shaspider"="Spider","Shacara"="Carabid","Shabird"="Bird","Shabat"="Bat"))



grid.col = c(rep("palegreen",8))
lty_fd <- cut(null_f$P, breaks = c(0,0.8,0.9,0.95,1),labels = c(0,3,2,1), include.lowest = TRUE)
lty_fd <-  as.numeric(levels(lty_fd))[lty_fd]
lty_lwd <- ifelse(lty_fd == 1, 3, ifelse(lty_fd == 2, 1, ifelse(lty_fd == 3, 0.5, 0.1)))


png("figures/observed_d.png",width=1000,height=800,pointsize = 18)
par(cex = 1.3)
circos.par(canvas.xlim=c(-1,1.6))
chordDiagram(null_d[,1:3], grid.col = grid.col, 
             annotationTrack = c("grid", "axis"), annotationTrackHeight = uh(9, "mm"),
             col=colorRamp2(seq(-0.25,0.25,length=8),brewer.pal(8,"BrBG")),
             link.lty = lty_fd, link.border = "black", link.lwd = lty_lwd)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black")
}

rect(1.1,0.2,1.2,0.9)
rasterImage(as.raster(matrix(brewer.pal(8,"BrBG")[8:1])),1.1,0.2,1.2,0.9)
text(x=1.25,y=1,labels="Correlation")
text(x=1.35,y=seq(0.2,0.9,length=4),labels = c("-0.25","-0.12","0.12","0.25"))
segments(rep(1.2,4),seq(0.2,0.9,length=4),rep(1.25,4),seq(0.2,0.9,length=4))
mtext("Observed correlation diagramm diversity - diversity")
text(1.4,0,labels = "P(|Correlation| > 0)")
segments(rep(1.1,4),seq(-0.4,-0.1,length=4),rep(1.2,4),seq(-0.4,-0.1,length=4),lty = c(0,3,2,1),lwd = c(0.1, 0.5, 1, 3))
text(x=rep(1.35,4),y=seq(-0.4,-0.1,length=4),labels = c("<0.80","0.80 - 0.90","0.90 - 0.95",">0.95"))
circos.clear()
dev.off()

# done!

# now look at what happens after model fitting

p_mult <- posterior_samples(m_mult)
p_mult[,161:280] %>%
  gather("corr","value") %>%
  group_by(corr) %>%
  summarise(R_mult = median(value), P_mult = sum(value > 0) / 4000) %>%
  mutate(P_mult = ifelse(P_mult < 0.5, 1 - P_mult, P_mult)) %>%
  arrange(desc(P_mult)) %>%
  separate(corr,c("drop","From","To")) %>%
  select(-drop) %>%
  mutate(Type_from = ifelse(From %in% c("Cstock", "Decomp", "Biomass", "Pgerm", "Herbivory", "Predation", "Birdsmi", "frass"),"function","diversity")) %>%
  mutate(Type_to = ifelse(To %in% c("Cstock", "Decomp", "Biomass", "Pgerm", "Herbivory", "Predation", "Birdsmi", "frass"),"function","diversity")) -> mult_dd

dd_all <- left_join(null_dd,mult_dd,by=c("From","To"))
dd_all$change <- "none"
dd_all$change <- ifelse(dd_all$P > 0.9 & dd_all$P_mult < 0.9, "extrinsic",dd_all$change)
dd_all$change <- ifelse(dd_all$P > 0.9 & dd_all$P_mult > 0.9, "intrinsic",dd_all$change)
dd_all$change <- ifelse(dd_all$P < 0.9 & dd_all$P_mult > 0.9, "appearing",dd_all$change)

# cool stuff now plot this on a chordDiagramm
dd_all$change_type <- paste(dd_all$Type_from.x, dd_all$Type_to.x,sep="_")
dd_all$change_type[dd_all$change_type == "diversity_function"] <- "function_diversity"
# first function - diversity
dd_all %>%
  filter(change_type == "function_diversity") %>%
  mutate(From = factor(From), To = factor(To)) %>%
  filter(change != "none") -> dd_fd

dd_fd$From <- revalue(dd_fd$From, c("Biomass"="TBiom","Birdsmi"="BBiom","frass"="Insect biomass"))
dd_fd$To <- revalue(dd_fd$To, c("Shaveg"="Vegetation","Shaherb"="Herbivore","Shaiso"="Isopod","Shadip"="Dipl.","Shaspider"="Spider","Shacara"="Carabid","Shabird"="Bird","Shabat"="Bat"))

link_col <- ifelse(dd_fd$change == "intrinsic","red",ifelse(dd_fd$change == "extrinsic","blue","green"))

grid.col = c(rep("salmon",7),rep("palegreen",6))

png("figures/changing_fd.png",width=1000,height=800,pointsize = 18)
par(cex = 1.3)
circos.par(canvas.xlim=c(-1,1.6))
chordDiagram(dd_fd[,c(1,2,3)], grid.col = grid.col, 
             annotationTrack = c("grid", "axis"), annotationTrackHeight = uh(9, "mm"),
             col=link_col, link.border = "black")
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black")
}

rect(1.1,0.2,1.2,0.3,col = "red")
rect(1.1,0.35,1.2,0.45,col = "green")
rect(1.1,0.5,1.2,0.6,col = "blue")
text(1.25, 0.25, labels = "Intrinsic", adj = 0.1)
text(1.25, 0.4, labels = "Appearing", adj = 0.1)
text(1.25, 0.55, labels = "Extrinsic", adj = 0.1)
mtext("Change in function - diversity correlation")
circos.clear()
dev.off()

# now function - function
dd_all %>%
  filter(change_type == "function_function") %>%
  mutate(From = factor(From), To = factor(To)) %>%
  filter(change != "none") -> dd_f

dd_f$From <- revalue(dd_f$From, c("Biomass"="TBiom","Birdsmi"="BBiom","frass"="Insect biomass"))
dd_f$To <- revalue(dd_f$To, c("Biomass"="TBiom","Birdsmi"="BBiom","frass"="Insect biomass"))

link_col <- ifelse(dd_f$change == "intrinsic","red",ifelse(dd_f$change == "extrinsic","blue","green"))

grid.col = c(rep("salmon",8))

png("figures/changing_f.png",width=1000,height=800,pointsize = 18)
par(cex = 1.3)
circos.par(canvas.xlim=c(-1,1.6))
chordDiagram(dd_f[,c(1,2,3)], grid.col = grid.col, 
             annotationTrack = c("grid", "axis"), annotationTrackHeight = uh(9, "mm"),
             col=link_col, link.border = "black")
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black")
}

rect(1.1,0.2,1.2,0.3,col = "red")
rect(1.1,0.35,1.2,0.45,col = "green")
rect(1.1,0.5,1.2,0.6,col = "blue")
text(1.25, 0.25, labels = "Intrinsic", adj = 0.1)
text(1.25, 0.4, labels = "Appearing", adj = 0.1)
text(1.25, 0.55, labels = "Extrinsic", adj = 0.1)
mtext("Change in function - function correlation")
circos.clear()
dev.off()

# now diversity - diversity
dd_all %>%
  filter(change_type == "diversity_diversity") %>%
  mutate(From = factor(From), To = factor(To)) %>%
  filter(change != "none") -> dd_d

dd_d$To <- revalue(dd_d$To, c("Shaveg"="Vegetation","Shaherb"="Herbivore","Shaiso"="Isopod","Shadip"="Dipl.","Shaspider"="Spider","Shacara"="Carabid","Shabird"="Bird","Shabat"="Bat"))
dd_d$From <- revalue(dd_d$From, c("Shaveg"="Vegetation","Shaherb"="Herbivore","Shaiso"="Isopod","Shadip"="Dipl.","Shaspider"="Spider","Shacara"="Carabid","Shabird"="Bird","Shabat"="Bat"))

link_col <- ifelse(dd_d$change == "intrinsic","red",ifelse(dd_d$change == "extrinsic","blue","green"))

grid.col = c(rep("palegreen",5))

png("figures/changing_d.png",width=1000,height=800,pointsize = 18)
par(cex = 1.3)
circos.par(canvas.xlim=c(-1,1.6))
chordDiagram(dd_d[,c(1,2,3)], grid.col = grid.col, 
             annotationTrack = c("grid", "axis"), annotationTrackHeight = uh(9, "mm"),
             col=link_col, link.border = "black")
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black")
}

rect(1.1,0.2,1.2,0.3,col = "red")
rect(1.1,0.35,1.2,0.45,col = "green")
rect(1.1,0.5,1.2,0.6,col = "blue")
text(1.25, 0.25, labels = "Intrinsic", adj = 0.1)
text(1.25, 0.4, labels = "Appearing", adj = 0.1)
text(1.25, 0.55, labels = "Extrinsic", adj = 0.1)
mtext("Change in diversity - diversity correlation")
circos.clear()
dev.off()



dd_all %>%
  group_by(change_type) %>%
  summarise(Ex_n = sum(change == "extrinsic"),
            In_n = sum(change == "intrinsic"),
            Ap_n = sum(change == "appearing")) %>%
  gather("direction","number",-change_type) %>%
  group_by(change_type) %>%
  mutate(prop = number / sum(number)) %>%
  arrange(change_type) -> dd_2

gg_cor <- ggplot(dd_2,aes(x=direction,y=prop,fill=change_type)) +
  geom_bar(stat="identity",position="dodge") +
  scale_x_discrete(labels = c("Appearing correlation","Extrinsic correlation","Intrinsic correlation"), name = "Type of correlation") +
  labs(y = "Proportion of correlation changes")

ggsave("figures/correlation_shift.png",gg_cor)

# 4. desirability

## helper functions
### this function reverse variables that need to be minimized
revv <- function(vec){
  return(- vec + (max(vec) + min(vec)))
}

desirab <- function(linpred, direction, importance, quantiles = c(0.1, 0.5, 0.9), sum = FALSE){
  for(i in 1:dim(linpred)[3]){ # reverse the functions that need to
    if(direction[i] == "minimize"){
      linpred[,,i] <- t(apply(linpred[,,i], 1, revv))
    }
  }
  
  avg_fw <- apply(linpred[,,1:8], c(1, 2), function(x) sum(importance[1:8] * x) / sum(importance[1:8]))
  if(sum){
    avg_fs <- rowSums(avg_fw)
    avg_fqq <- quantile(avg_fs, probs = quantiles)
    avg_fout <- data.frame(LCI = avg_fqq[1], Median = avg_fqq[2], UCI = avg_fqq[3], type = "functioning")
  }
  else{
    avg_fqq <- apply(avg_fw, 2, quantile, probs = quantiles) 
    avg_fout <- data.frame(LCI = avg_fqq[1,], Median = avg_fqq[2,], UCI = avg_fqq[3,], type = "functioning")
  }
  
  
  avg_dw <- apply(linpred[,,9:16], c(1, 2), function(x) sum(importance[9:16] * x) / sum(importance[9:16]))
  if(sum){
    avg_ds <- rowSums(avg_dw)
    avg_dqq <- quantile(avg_ds, probs = quantiles)
    avg_dout <- data.frame(LCI = avg_dqq[1], Median = avg_dqq[2], UCI = avg_dqq[3], type = "diversity")
  }
  else{
    avg_dqq <- apply(avg_dw, 2, quantile, probs = quantiles) 
    avg_dout <- data.frame(LCI = avg_dqq[1,], Median = avg_dqq[2,], UCI = avg_dqq[3,], type = "diversity")
  }
  avg_out <- rbind(avg_fout, avg_dout)
  
  
  return(avg_out)
}

## 4.0 import importance scores
imp <- read.csv("data/importance_score.csv")

## 4.1 desirability at plot-scale
newdat_plot <- expand.grid(speccomb = c("fsyl","qrob","qrub","all"), prox.std = c(0, -1.37), edge.std = c(2, -1.1), dens.std = 0)[c(1:4,13:16),]
newdat_plot <- expand.grid(speccomb = c("fsyl","qrob","qrub","all"), prox.std = 0, edge.std = c(2, -1.1), dens.std = 0)
newdat_plot <- expand.grid(speccomb = c("fsyl","qrob","qrub","all"), prox.std = c(-1.37, 3.5), edge.std = 0, dens.std = 0)


# get predictions
pred_plot <- posterior_linpred(m_mult, newdata = newdat_plot)
# get desirability scores
dd_plotf <- desirab(pred_plot, imp$direction_f, imp$weight_f)
dd_plotf$perspective <- "productivist"

dd_plotc <- desirab(pred_plot, imp$direction_c, imp$weight_c)
dd_plotc$perspective <- "conservationist"

# dd_plota <- desirab(pred_plot, rep("a",16), rep(5, 16))
# dd_plota$perspective <- "average"

dd_plot <- rbind(dd_plotf, dd_plotc)
dd_plot$fragmentation <- rep(rep(c("high frag.","low frag."),each = 4), 4)
dd_plot$speccomb <- rep(newdat_plot$speccomb, 4)
#dd_plot$X <- rep(1:4, 8) + rep(c(-0.1, 0, 0.1), each = 8)
dd_plot$X2 <- rep(1:4, 8) + c(rep(rep(c(-0.15, 0.05), each = 4), 2),rep(rep(c(-0.05, 0.15), each = 4), 2))

gg_plot <- ggplot(dd_plot, aes(x = X2, y = Median, ymin = LCI, ymax = UCI, color = fragmentation, shape = perspective)) +
  geom_linerange() +
  geom_point(size = 2.5) +
  facet_grid(type ~ .) +
  scale_x_continuous(breaks = 1:4, labels = c("fsyl","qrob","qrub","all")) +
  labs(x = "Tree composition", y = "Desirability index (with 80% CrI)")
  
ggsave("figures/desirab_plot.png",gg_plot)

dd_plot$top <- "plot-level"
# to combine with landscape-scale
gg_plot <- ggplot(dd_plot, aes(x = X2, y = Median, ymin = LCI, ymax = UCI, color = fragmentation, shape = perspective)) +
  geom_linerange() +
  geom_point(size = 2.5) +
  facet_grid(type ~ top) +
  scale_x_continuous(breaks = 1:4, labels = c("fsyl","qrob","qrub","all")) +
  labs(x = "", y = "Desirability scores (with 80% CrI)") +
  theme(legend.position = "none", strip.text.y = element_blank()) 
  

# desirability at the landscape scale
## forester
newdat_land1 <- data.frame(speccomb = rep(c("fsyl","qrob","qrub"), times = c(17,18,18)), edge.std = 2, prox.std = 0, dens.std = 0)
newdat_land2 <- data.frame(speccomb = rep("all", 53), edge.std = 2, prox.std = 0, dens.std = 0)
newdat_land3 <- data.frame(speccomb = rep(c("fsyl","qrob","qrub"), times = c(17,18,18)), edge.std = -1.11, prox.std = -1.37, dens.std = 0)
newdat_land4 <- data.frame(speccomb = rep("all", 53), edge.std = -1.11, prox.std = -1.37, dens.std = 0)

pred_1 <- posterior_linpred(m_mult, newdata = newdat_land1)
pred_2 <- posterior_linpred(m_mult, newdata = newdat_land2)
pred_3 <- posterior_linpred(m_mult, newdata = newdat_land3)
pred_4 <- posterior_linpred(m_mult, newdata = newdat_land4)
pred_l <- list("monoculture\nhigh frag."=pred_1, "mixtures\nhigh frag." = pred_2, "monoculture\nlow frag." = pred_3, "mixtures\nlow frag." = pred_4)


dd_f <- ldply(pred_l, function(x) desirab(x, imp$direction_f, imp$weight_f, sum = TRUE))
dd_c <- ldply(pred_l, function(x) desirab(x, imp$direction_c, imp$weight_c, sum = TRUE))
dd_f$perspective <- "productivist"
dd_c$perspective <- "conservationist"

## observed
pred_o <- posterior_linpred(m_mult)
dd_o <- desirab(pred_o, imp$direction_f, imp$weight_f, sum = TRUE)
dd_o$.id <- "original"
dd_o$perspective <- "forester"
dd_o <- dd_o[,c(5,1:4,6)]
dd_o <- rbind(dd_o,dd_o)
dd_o$perspective[3:4] <- "conservationist"

dd_all <- rbind(dd_f,dd_c)

dd_all <- separate(dd_all, ".id", c("tree_composition", "fragmentation"), sep = "\n")
dd_all$X <- rep(1:2, each = 2) + rep(c(-0.15, 0.05, -0.05,0.15), each = 4)
dd_all$truc <- "landscape"

gg_land <- ggplot(dd_all,aes(x=X,y=Median, ymin = LCI, ymax = UCI, color = fragmentation, shape = perspective)) +
  geom_linerange() +
  geom_point(size = 2.5) +
  facet_grid(type ~ .) +
  labs(x = "Scenario (tree composition / fragmentation)", y = "Desirability index (with 80% CrI)") +
  scale_x_continuous(breaks = 1:2, labels = c("monoculture","mixture"), limits = c(0.5,2.5))

ggsave("figures/desirab_landscape.png", gg_land)

# together with plot scale
gg_land <- ggplot(dd_all,aes(x=X,y=Median, ymin = LCI, ymax = UCI, color = fragmentation, shape = perspective)) +
  geom_linerange() +
  geom_point(size = 2.5) +
  facet_grid(type ~ truc) +
  labs(x = "", y = "") +
  scale_x_continuous(breaks = 1:2, labels = c("monoculture","mixture"), limits = c(0.5,2.5))

gg_desirab <- grid.arrange(gg_plot, gg_land, bottom = "Tree species composition", widths = c(6.5, 8.5))

ggsave("figures/03_desirab.png", gg_desirab)

# 5. Radar plot, another way to show trade-offs and synergies
newdat <- data.frame(speccomb = c("fsyl", "qrob", "qrub", "all"), edge.std = 0, prox.std = 0, dens.std = 0)

pred_comp <- posterior_linpred(m_mult, newdata = newdat)
pred_comp <- adply(pred_comp,c(2,3),quantile, probs=c(0.1,0.5,0.9))
pred_comp$speccomb <- rep(c("fsyl","qrob","qrub","all"), 16)
names(pred_comp)[3:5] <- c("LCI","Median","UCI")
tmp <- melt(pred_comp, id.vars = c(2,6),measure.vars = 3:5)

# now for function
pred_r <- dcast(tmp[c(1:32,65:96,129:160),],speccomb + variable ~ X2, value.var = "value")
rownames(pred_r) <- paste(pred_r$speccomb,pred_r$variable)
pred_r <- pred_r[,-c(1,2)]
pred_r <- rbind(rep(1.75, 8), rep(-1.75,8),rep(0,8), pred_r)
colnames(pred_r) <- c("C stock", "Decomposition", "Tree\nbiom.", "Regeneration", "Herbivory", "Predation", "Bird\nbiom.", "Insect biomass")
pred_r <- pred_r[,c(1:4, 8, 5, 7, 6)]



# now for diversity
pred_d <- dcast(tmp[c(33:64,97:128,161:192),],speccomb + variable ~ X2, value.var = "value")
rownames(pred_d) <- paste(pred_d$speccomb,pred_d$variable)
pred_d <- pred_d[,-c(1,2)]
pred_d <- rbind(rep(1.75, 8), rep(-1.75,8),rep(0,8), pred_d)
colnames(pred_d) <- c("Vegetation", "Herbivore", "Cara.", "Spider", "Isopod", "Diplopod", "Bird", "Bat")

cols_border <- c("black",rep(viridis(4),each=3))
cols_in <- NA


# all in one
png("figures/spid_comp.png", width = 1200, height = 800)
m <- matrix(c(1, 2, 3, 3), ncol = 2, byrow = TRUE)
layout(m, heights = c(0.8, 0.2))

par(mar = c(1,1,1,1), cex = 1.5)
radarchart(pred_r, axistype = 1, pcol = cols_border, pfcol = cols_in, plwd = c(4,rep(c(2,4,2),8)), plty = c(1,rep(c(2,1,2),8)),pty = c(NA,rep(c(NA,16,NA),8)),
           cglcol = "grey", cglty = 1, axislabcol = "black", caxislabels = seq(-1.5, 1.5, length = 5), cglwd = 0.8,vlcex = 1,
           title = "Ecosystem function")

radarchart(pred_d, axistype = 1, pcol = cols_border, pfcol = cols_in, plwd = c(4,rep(c(2,4,2),8)), plty = c(1,rep(c(2,1,2),8)),pty = c(NA,rep(c(NA,16,NA),8)),
           cglcol = "grey", cglty = 1, axislabcol = "black", caxislabels = seq(-1.5, 1.5, length = 5), cglwd = 0.8, vlcex = 1,
           title = "Diversity")
par(mar = c(0,0,0,0))

plot(1, type = "n", axes = FALSE, xlab="",ylab="")


legend("top", legend = c("all", "fsyl", "qrob", "qrub"), bty = "n", col = viridis(4), pch = 16, lty = 1, lwd = 4, horiz = TRUE,
       title = "Tree composition", cex = 1.5)

dev.off()

# now for edge effects
newdat <- expand.grid(speccomb = c("fsyl","qrob","qrub","all"), edge.std = c(-1.2, 2), prox.std = 0, dens.std = 0)

pred_e <- posterior_linpred(m_mult, newdata = newdat)
pred_e2 <- adply(pred_e,c(1,3), function(x) data.frame(Low = mean(x[1:4]), High = mean(x[5:8])))
pred_e2 %>%
  group_by(X2) %>%
  summarise(LowM = median(Low), HighM = median(High), LowLCI = quantile(Low, probs = 0.1), HighLCI = quantile(High,probs = 0.1), 
            LowUCI = quantile(Low, probs = 0.9), HighUCI = quantile(High,probs = 0.9)) -> pred_e3

pred_f <- t(pred_e3[1:8,2:7])
colnames(pred_f) <- c("C stock", "Decomposition", "Tree\nbiom.", "Regeneration","Herbivory","Predation","Bird\nbiom.", "Insect biomass")
pred_f <- as.data.frame(rbind(rep(1.25,8),rep(-1.25,8),rep(0,8),pred_f))
pred_f <- pred_f[,c(1:4, 8, 5, 7, 6)]

pred_d <- t(pred_e3[9:16,2:7])
colnames(pred_d) <- c("Vegetation", "Herbivore", "Cara.", "Spider","Isopod","Diplopod","Bird", "Bat")
pred_d <- as.data.frame(rbind(rep(1.25,8),rep(-1.25,8),rep(0,8),pred_d))

cols_bordere <- c("black",rep(c("orange","royalblue"),3))
cols_ine <- NA

png("figures/spid_edge.png", width = 1200, height = 800)
m <- matrix(c(1, 2, 3, 3), ncol = 2, byrow = TRUE)
layout(m, heights = c(0.8, 0.2))

par(mar = c(1,1,1,1), cex = 1.5)

radarchart(pred_f, axistype = 1, pcol = cols_bordere, pfcol = cols_ine, plwd = c(4,4,4,2,2,2,2), plty = c(1,1,1,2,2,2,2), pty = c(NA,16,16,NA,NA,NA,NA),
           cglcol = "grey", cglty = 1, axislabcol = "black", caxislabels = seq(-1.25, 1.25, length = 5), cglwd = 1,vlcex = 1,
           title = "")

radarchart(pred_d, axistype = 1, pcol = cols_bordere, pfcol = cols_ine, plwd = c(4,4,4,2,2,2,2), plty = c(1,1,1,2,2,2,2),pty = c(NA,16,16,NA,NA,NA,NA),
           cglcol = "grey", cglty = 1, axislabcol = "black", caxislabels = seq(-1.25, 1.25, length = 5), cglwd = 1,vlcex = 1,
           title = "")

par(mar = c(0,0,0,0))

plot(1, type = "n", axes = FALSE, xlab="",ylab="")


legend("top", legend = c("Low", "High"), bty = "n", col = cols_bordere[2:3], pch = 16, lty = 1, lwd = 4, horiz = TRUE,
       title = "Amount of edges:",cex=1.5)

dev.off()

# now for proximity
newdat <- expand.grid(speccomb = c("fsyl","qrob","qrub","all"), edge.std = 0, prox.std = c(-1.5, 3.5), dens.std = 0)

pred_p <- posterior_linpred(m_mult, newdata = newdat)
pred_p2 <- adply(pred_p,c(1,3), function(x) data.frame(Low = mean(x[1:4]), High = mean(x[5:8])))
pred_p2 %>%
  group_by(X2) %>%
  summarise(LowM = median(Low), HighM = median(High), LowLCI = quantile(Low, probs = 0.1), HighLCI = quantile(High,probs = 0.1), 
            LowUCI = quantile(Low, probs = 0.9), HighUCI = quantile(High,probs = 0.9)) -> pred_p3

pred_f <- t(pred_p3[1:8,2:7])
colnames(pred_f) <- c("C stock", "Decomposition", "Tree\nbiom.", "Regeneration","Herbivory","Predation","Bird\nbiom.", "Insect biomass")
pred_f <- as.data.frame(rbind(rep(3,8),rep(-3,8),rep(0,8),pred_f))
pred_f <- pred_f[,c(1:4, 8, 5, 7, 6)]

pred_d <- t(pred_p3[9:16,2:7])
colnames(pred_d) <- c("Vegetation", "Herbivore", "Cara.", "Spider","Isopod","Diplopod","Bird", "Bat")
pred_d <- as.data.frame(rbind(rep(3,8),rep(-3,8),rep(0,8),pred_d))

cols_borderp <- c("black",rep(c("peru","slateblue"),3))
cols_inp <- NA

png("figures/spid_prox.png", width = 1200, height = 800)
m <- matrix(c(1, 2, 3, 3), ncol = 2, byrow = TRUE)
layout(m, heights = c(0.8, 0.2))

par(mar = c(1,1,1,1), cex = 1.5)

radarchart(pred_f, axistype = 1, pcol = cols_borderp, pfcol = cols_inp, plwd = c(4,4,4,2,2,2,2), plty = c(1,1,1,2,2,2,2), pty = c(NA,16,16,NA,NA,NA,NA),
           cglcol = "grey", cglty = 1, axislabcol = "black", caxislabels = seq(-3, 3, length = 5), cglwd = 1,vlcex = 1,
           title = "")

radarchart(pred_d, axistype = 1, pcol = cols_borderp, pfcol = cols_inp, plwd = c(4,4,4,2,2,2,2), plty = c(1,1,1,2,2,2,2),pty = c(NA,16,16,NA,NA,NA,NA),
           cglcol = "grey", cglty = 1, axislabcol = "black", caxislabels = seq(-3, 3, length = 5), cglwd = 1,vlcex = 1,
           title = "")

par(mar = c(0,0,0,0))

plot(1, type = "n", axes = FALSE, xlab="",ylab="")


legend("top", legend = c("Low", "High"), bty = "n", col = cols_borderp[2:3], pch = 16, lty = 1, lwd = 4, horiz = TRUE,
       title = "Proximity to other forest fragments:",cex=1.5)

dev.off()


### unrelated stuff for the soiltemp meeting
dd <- cbind(func[,c("Decomp","Herbivory","Predation")], pred_dat[,c(3,4)])
names(dd)[4] <- "edge"
dum <- read.csv("~/Documents/dummy.csv")
dum %>%
group_by(id_plot) %>%
summarise(N = n(), S = sum(attacked, na.rm=TRUE))->dum_dd
dd$N <- dum_dd$N
dd$S <- dum_dd$S


m_d <- lm(Decomp ~ edge * speccomb, dd)
newdat <- expand.grid(edge = c(0, 500), speccomb = c("all", "fsyl", "qrob","qrub", "fsyl_qrob", "fsyl_qrub", "qrob_qrub"))
pedd <- predict(m_d, newdata = newdat, se.fit = TRUE)
newdat$Decomp <- pedd$fit
newdat$LCI <- newdat$Decomp - 1.96 * pedd$se.fit
newdat$UCI <- newdat$Decomp + 1.96 * pedd$se.fit

ggplot(dd, aes(x = edge, y = Decomp, color = speccomb, fill = speccomb)) +
  geom_point() +
  geom_line(data = newdat) +
  geom_ribbon(data = newdat, aes(ymin = LCI, ymax = UCI), alpha = 0.1) +
  facet_wrap(~speccomb) +
  labs(x = "Edge length in a buffer of 100m", y= "Decomposition rate") +
  scale_fill_discrete(name = "Tree species\ncomposition") +
  scale_color_discrete(name = "Tree species\ncomposition") +
  theme(text = element_text(size=20))

# now for herbivory
m_h <- lm(Herbivory ~ edge * speccomb, dd)

ggplot(dd, aes(x = edge, y = Herbivory, color = speccomb, fill = speccomb)) +
  geom_point() +
  facet_wrap(~speccomb) +
  labs(x = "Edge length in a buffer of 100m", y= "Herbivory rate") +
  scale_fill_discrete(name = "Tree species\ncomposition") +
  scale_color_discrete(name = "Tree species\ncomposition") +
  theme(text = element_text(size=20))

# now for predation
m_p <- glm(cbind(S, N - S) ~ edge * speccomb, dd, family = "binomial")

newdat <- expand.grid(edge = seq(0,500, length = 20), speccomb = c("all", "fsyl", "qrob","qrub", "fsyl_qrob", "fsyl_qrub", "qrob_qrub"))
pedd <- predict(m_p, newdata = newdat, se.fit = TRUE, type = "response")
newdat$Predation <- pedd$fit
newdat$LCI <- newdat$Predation - 1.96 * pedd$se.fit
newdat$UCI <- newdat$Predation + 1.96 * pedd$se.fit

ggplot(dd, aes(x = edge, y = Predation, color = speccomb, fill = speccomb)) +
  geom_point() +
  geom_line(data = newdat) +
  geom_ribbon(data = newdat, aes(ymin = LCI, ymax = UCI), alpha = 0.1) +
  facet_wrap(~speccomb) +
  labs(x = "Edge length in a buffer of 100m", y= "Predation rate") +
  scale_fill_discrete(name = "Tree species\ncomposition") +
  scale_color_discrete(name = "Tree species\ncomposition")+
  theme(text = element_text(size=20))

# coming back to variance explained
pp <- posterior_samples(m_mult)

# dev, focus on one variable at a time
pp2 <- pp[,c(1,17:24)] # all coefs for C stock



# make a function
get_sd <- function(post, name = "Cstock"){
  # grab fixed effects
  fixx <- post[,grep(paste0("b_",name), names(post))]
  # the sd for speccomb
  comp_sd <- apply(fixx[,1:7],1,function(b) sd(b[as.numeric(brm_dat$speccomb)]))
  # the sd for edge
  edge_sd <- abs(fixx[,9])
  # the sd for prox
  prox_sd <- abs(fixx[,8])
  # the sd for density
  # dens_sd <- abs(fixx[,10])
  # the sd of the resid
  res_sd <- resid_all[,name]
  sd.all <- cbind(comp_sd, edge_sd, prox_sd, res_sd)
  out <- tidyMCMC(100 * sd.all / rowSums(sd.all), conf.int = TRUE, conf.level = 0.8)
  out$variable <- name
  out$fraction <- c("composition","edge","proximity", "residual")
  return(out)
}

# go throught all the functions
nn <- names(brm_dat)[2:17]
nn <- gsub("\\.","",nn)
nn <- nn[c(1:4,8,5,7,6,9:16)]
post <- posterior_samples(m_mult)
residd <- residuals(m_mult, summary = FALSE)
resid_all <- apply(residd, c(1,3), sd)

sd_all <- ldply(nn, function(n) get_sd(post, n))
sd_all$X <- rep(1:16, each = 4) + rep(c(-0.15, -0.07, 0.07, 0.15), 16)

sd_all %>%
  group_by(fraction) %>%
  summarise(estimate = mean(estimate)) -> sd_avg

lbls <- c("C. stock", "Decomposition", "Tree biomass", "Regeneration", "Insect biomass", "Herbivory", "Bird biomass", "Predation",
          "Vegetation div.", "Herbivore div.", "Carabid div.", "Spider div.", "Isopod div.", "Diplopod div.", "Bird div.", "Bat div.")

gg_sd <- ggplot(sd_all,aes(x= X, y = estimate, ymin = conf.low, ymax = conf.high, color = fraction)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  scale_x_continuous(breaks = 1:16, labels = lbls) +
  geom_hline(data = sd_avg, aes(yintercept = estimate, color = fraction), linetype = "dashed") +
  labs(x="", y = "Percentage of explained variance") +
  scale_color_discrete(name = "Effects")

ggsave("figures/explained.png", gg_sd, width = 8, height = 7)

# the sd for speccomb
comp_sd <- apply(pp2[,1:7],1,function(b) sd(b[as.numeric(pred_dat$speccomb)]))
# the sd for edge
edge_sd <- abs(pp2$b_Cstock_edge.std)
# the sd for prox
prox_sd <- abs(pp2$b_Cstock_prox.std)
# the sd of the resid
res_sd <- apply(residd[,,1],1,sd)
sd.all <- cbind(comp_sd, edge_sd, prox_sd, res_sd)
out <- tidyMCMC(100 * sd.all / rowSums(sd.all), conf.int = TRUE, conf.level = 0.8)

# TX Make summary table of dataset
rawdat <- cbind(func[,c("C_stock", "Decomp","Biomass", "P_germ", "frass", "Herbivory", "Bird_smi", "Predation")], div_all[,seq(3, 33, 4)])
rawdat$Decomp <- rawdat$Decomp * 100
rawdat$Predation <- rawdat$Predation * 100
rawdat$P_germ <- rawdat$P_germ * 100
rawdat$frass <- rawdat$frass / 2.25
tt <- t(apply(rawdat, 2, function(x) c(paste0(round(mean(x, na.rm = TRUE), 2), " (", round(sd(x, na.rm = TRUE), 2), ")"),
                               paste0(round(min(x, na.rm = TRUE), 2), " - ", round(max(x, na.rm = TRUE), 0)))))

xtable(tt)

# figures for the presentation
imp <- read.csv("data/importance_score.csv")

ff <- factor(c("Carbon", "Decomposition", "Tree biom.", "Regeneration", "Herbivory", "Predation", "Bird biom.", "Insect biom.",
               "Vegetation", "Herbivore", "Carabid", "Spider", "Isopod", "Diplopod", "Bird", "Bat"),
             levels = c("Carbon", "Decomposition", "Tree biom.", "Regeneration", "Insect biom.", "Herbivory", "Bird biom.", "Predation",
                        "Vegetation", "Herbivore", "Carabid", "Spider", "Isopod", "Diplopod", "Bird", "Bat"))

imp %>%
  rename(productivist = weight_f, conservationist = weight_c, variable = function.) %>%
  gather("perspective", "weight", productivist, conservationist) %>%
  mutate(variable = rep(ff, 2)) -> dd

ggplot(dd, aes(x = variable, y = weight, fill = perspective)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust = 0.8)) +
  geom_vline(xintercept = 8.5, linetype = "dashed", width = 1.5) +
  labs(x="")
  
# some formatting of carabids data
cara %>%
  separate("CODE.Pallieter", c("Plot", "Date", "Pitfall_id"), sep = "-") %>%
  rename(Code = Species) %>%
  left_join(sp_code, by = "Code") %>%
  select(Plot, Date, Pitfall_id, Code, Genus, species, SEX, Aantal) %>%
  arrange(Plot, Date, Pitfall_id, Code) %>%
  filter(!is.na(Date)) %>%
  filter(species != "") -> cara_dd

write.table(cara_dd,"treeweb_carabid_format.csv", sep = ",", row.names = FALSE)


pred_ex <- pred_f[1:4,]
pred_ex[4,] <- c(3, 1, 0, 0, -1, -3, -2, 2)

cols_ex <- c("black", "red")

radarchart(pred_ex, axistype = 1, pcol = cols_ex, pfcol = cols_inp, plwd = c(4,4), plty = c(1,1),
                    cglcol = "grey", cglty = 1, axislabcol = "black", caxislabels = seq(-3, 3, length = 5), cglwd = 1,vlcex = 1,
                    title = "", pty = c(NA, 16))



# try again at another time ...
con <- DBI::dbConnect(RPostgreSQL::PostgreSQL(),
dbname = "gis",
user = "hertzog",
host = "134.110.32.102",
port = "14001",
password = rstudioapi::askForPassword("Database password")
)

